use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Seek, SeekFrom};
use byteorder::{BigEndian, LittleEndian, ReadBytesExt};
use minifb::{Key, MouseButton, MouseMode, Window, WindowOptions};

const WIDTH: usize = 800;
const HEIGHT: usize = 600;

#[derive(Debug)]
enum Shape {
    Point(f64, f64),
    Poly(Vec<(f64, f64)>),
}

// نرمالایز کردن مختصات جغرافیایی به مختصات صفحه نمایش
fn normalize(x: f64, y: f64, view: (f64, f64, f64, f64)) -> (usize, usize) {
    let (xmin, xmax, ymin, ymax) = view;
    let px = (((x - xmin) / (xmax - xmin)) * (WIDTH as f64)) as usize;
    let py = (((1.0 - (y - ymin) / (ymax - ymin))) * (HEIGHT as f64)) as usize;
    (px.min(WIDTH - 1), py.min(HEIGHT - 1))
}

fn draw_point(buf: &mut [u32], x: usize, y: usize, color: u32) {
    if x < WIDTH && y < HEIGHT {
        buf[y * WIDTH + x] = color;
    }
}

// رسم خط بین دو نقطه (Bresenham)
fn draw_line(buf: &mut [u32], x0: usize, y0: usize, x1: usize, y1: usize, color: u32) {
    let dx = (x1 as isize - x0 as isize).abs();
    let dy = -(y1 as isize - y0 as isize).abs();
    let sx = if x0 < x1 { 1 } else { -1 };
    let sy = if y0 < y1 { 1 } else { -1 };
    let mut err = dx + dy;
    let (mut x, mut y) = (x0 as isize, y0 as isize);

    loop {
        if x >= 0 && y >= 0 && (x as usize) < WIDTH && (y as usize) < HEIGHT {
            buf[y as usize * WIDTH + x as usize] = color;
        }
        if x == x1 as isize && y == y1 as isize { break; }
        let e2 = 2 * err;
        if e2 >= dy { err += dy; x += sx; }
        if e2 <= dx { err += dx; y += sy; }
    }
}

// خواندن هندسه‌ها از فایل shp
fn read_shapes(path: &str) -> Result<(Vec<Shape>, (f64, f64, f64, f64)), Box<dyn std::error::Error>> {
    let mut file = BufReader::new(File::open(path)?);

    // خواندن باکس مرجع از هدر
    file.seek(SeekFrom::Start(36))?;
    let xmin = file.read_f64::<LittleEndian>()?;
    let ymin = file.read_f64::<LittleEndian>()?;
    let xmax = file.read_f64::<LittleEndian>()?;
    let ymax = file.read_f64::<LittleEndian>()?;
    file.seek(SeekFrom::Start(100))?;

    let mut shapes = Vec::new();

    loop {
        let pos = file.stream_position()?;
        if pos >= file.get_ref().metadata()?.len() {
            break;
        }
        // هر رکورد با شماره و طول محتوا شروع می‌شود
        let res = file.read_i32::<BigEndian>();
        if res.is_err() { break; }
        let _record_number = res?;
        let _content_length = file.read_i32::<BigEndian>()?;

        let shape_type = file.read_i32::<LittleEndian>()?;
        match shape_type {
            1 => { // Point
                let x = file.read_f64::<LittleEndian>()?;
                let y = file.read_f64::<LittleEndian>()?;
                shapes.push(Shape::Point(x, y));
            }
            3 => { // Polyline یا Polygon (هر دو مشابه)
                file.seek(SeekFrom::Current(32))?; // Bounding Box
                let num_parts = file.read_i32::<LittleEndian>()?;
                let num_points = file.read_i32::<LittleEndian>()?;

                let mut parts = Vec::new();
                for _ in 0..num_parts {
                    parts.push(file.read_i32::<LittleEndian>()?);
                }

                let mut points = Vec::new();
                for _ in 0..num_points {
                    let x = file.read_f64::<LittleEndian>()?;
                    let y = file.read_f64::<LittleEndian>()?;
                    points.push((x, y));
                }

                shapes.push(Shape::Poly(points));
            }
            _ => {
                // سایر انواع هندسه پشتیبانی نمی‌شود
                break;
            }
        }
    }

    Ok((shapes, (xmin, xmax, ymin, ymax)))
}

// خواندن رکوردهای DBF
fn read_dbf(path: &str) -> Result<Vec<HashMap<String, String>>, Box<dyn std::error::Error>> {
    let mut file = BufReader::new(File::open(path)?);

    // Read the 32-byte header
    let mut header = [0u8; 32];
    file.read_exact(&mut header)?;
    let num_records = u32::from_le_bytes([header[4], header[5], header[6], header[7]]) as usize;
    let header_len = u16::from_le_bytes([header[8], header[9]]) as usize;
    let record_len = u16::from_le_bytes([header[10], header[11]]) as usize;

    // Read field descriptors (32 bytes each until 0x0D)
    let mut fields = Vec::new();
    while file.stream_position()? < (header_len as u64) - 1 {
        let mut field = [0u8; 32];
        file.read_exact(&mut field)?;
        if field[0] == 0x0D {
            break;
        }
        let name = String::from_utf8_lossy(&field[0..11])
            .trim_end_matches('\0')
            .to_string();
        let length = field[16] as usize;
        fields.push((name, length));
    }

    // Skip field descriptor terminator
    file.seek(SeekFrom::Start(header_len as u64))?;

    // Read records
    let mut records = Vec::new();
    let mut buffer = vec![0u8; record_len];

    for _ in 0..num_records {
        match file.read_exact(&mut buffer) {
            Ok(_) => {
                if buffer[0] == 0x2A {
                    continue; // Deleted record
                }
                let mut map = HashMap::new();
                let mut offset = 1;
                for (name, len) in &fields {
                    let end = offset + *len;
                    if end > buffer.len() {
                        break;
                    }
                    let val = String::from_utf8_lossy(&buffer[offset..end])
                        .trim()
                        .to_string();
                    map.insert(name.clone(), val);
                    offset = end;
                }
                records.push(map);
            }
            Err(e) => {
                if e.kind() == std::io::ErrorKind::UnexpectedEof {
                    break; // پایان واقعی فایل
                } else {
                    return Err(Box::new(e));
                }
            }
        }
    }

    Ok(records)
}

// fn main() -> Result<(), Box<dyn std::error::Error>> {
//     let (shapes, bbox) = read_shapes("data.shp")?;
//     let dbf_records = read_dbf("data.dbf")?;

//     let mut window = Window::new("Shapefile Viewer", WIDTH, HEIGHT, WindowOptions::default())?;
//     let mut buffer = vec![0xFFFFFFu32; WIDTH * HEIGHT];

//     let (xmin, xmax, ymin, ymax) = bbox;

//     let mut selected_shape_idx: Option<usize> = None;

//     while window.is_open() && !window.is_key_down(Key::Escape) {
//         if let Some((mx, my)) = window.get_mouse_pos(MouseMode::Pass) {
//             if window.get_mouse_down(MouseButton::Left) {
//                 // تبدیل مختصات کلیک شده به مختصات shp
//                 let wx = xmin + (mx as f64 / WIDTH as f64) * (xmax - xmin);
//                 let wy = ymin + (1.0 - (my as f64 / HEIGHT as f64)) * (ymax - ymin);

//                 let mut closest_idx = None;
//                 let mut closest_dist = f64::MAX;

//                 for (i, shape) in shapes.iter().enumerate() {
//                     let dist = match shape {
//                         Shape::Point(x, y) => ((x - wx).powi(2) + (y - wy).powi(2)).sqrt(),
//                         Shape::Poly(pts) => pts.iter()
//                             .map(|(x, y)| ((x - wx).powi(2) + (y - wy).powi(2)).sqrt())
//                             .fold(f64::MAX, |a, b| a.min(b)),
//                     };
//                     if dist < closest_dist && dist < 0.01 * (xmax - xmin) {
//                         closest_dist = dist;
//                         closest_idx = Some(i);
//                     }
//                 }

//                 if let Some(idx) = closest_idx {
//                     selected_shape_idx = Some(idx);
//                     println!("Selected record index: {}", idx);
//                     if let Some(rec) = dbf_records.get(idx) {
//                         println!("DBF info:");
//                         for (k, v) in rec.iter() {
//                             println!("  {}: {}", k, v);
//                         }
//                     }
//                 } else {
//                     selected_shape_idx = None;
//                 }
//             }
//         }

//         // بک‌گراند سفید
//         buffer.fill(0xFFFFFF);

//         // رسم هندسه‌ها
//         for (i, shape) in shapes.iter().enumerate() {
//             let color = if Some(i) == selected_shape_idx { 0xFF0000 } else { 0x0000FF };
//             match shape {
//                 Shape::Point(x, y) => {
//                     let (px, py) = normalize(*x, *y, (xmin, xmax, ymin, ymax));
//                     draw_point(&mut buffer, px, py, color);
//                 }
//                 Shape::Poly(pts) => {
//                     for w in pts.windows(2) {
//                         let (x0, y0) = normalize(w[0].0, w[0].1, (xmin, xmax, ymin, ymax));
//                         let (x1, y1) = normalize(w[1].0, w[1].1, (xmin, xmax, ymin, ymax));
//                         draw_line(&mut buffer, x0, y0, x1, y1, color);
//                     }
//                 }
//             }
//         }

//         window.update_with_buffer(&buffer, WIDTH, HEIGHT)?;
//     }

//     Ok(())
// }

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let (shapes, _bbox_orig) = read_shapes("data.shp")?;
    let dbf_records = read_dbf("data.dbf")?;

    let mut window = Window::new("Shapefile Viewer", WIDTH, HEIGHT, WindowOptions::default())?;
    let mut buffer = vec![0xFFFFFFu32; WIDTH * HEIGHT];

    let mut view_bbox = _bbox_orig;

    let mut selected_shape_idx: Option<usize> = None;
    let mut mode = "identify"; // یا zoom

    // متغیرهای وضعیت موس برای pan
    let mut is_panning = false;
    let mut last_mouse_pos: Option<(f32, f32)> = None;

    println!("🎯 Mode: {}", mode);
    println!("🔁 Press 'Z' to toggle between Identify and Zoom");
    println!("➡️ Hold right mouse button and drag to Pan");

    while window.is_open() && !window.is_key_down(Key::Escape) {
        // تعویض بین حالت‌ها
        if window.is_key_pressed(Key::Z, minifb::KeyRepeat::No) {
            mode = if mode == "identify" { "zoom" } else { "identify" };
            println!("🎯 Mode changed to: {}", mode);
        }

        // pan با کلیک و درگ موس دکمه راست
        if window.get_mouse_down(MouseButton::Right) {
            if !is_panning {
                is_panning = true;
                last_mouse_pos = window.get_mouse_pos(MouseMode::Pass);
            } else if let Some((last_x, last_y)) = last_mouse_pos {
                if let Some((mx, my)) = window.get_mouse_pos(MouseMode::Pass) {
                    let dx = mx - last_x;
                    let dy = my - last_y;

                    let (_xmin, _xmax, _ymin, _ymax) = view_bbox;
                    let width = _xmax - _xmin;
                    let height = _ymax - _ymin;

                    // جابجایی view_bbox متناسب با حرکت موس
                    let dx_map = -dx as f64 / WIDTH as f64 * width;
                    let dy_map = dy as f64 / HEIGHT as f64 * height;

                    view_bbox = (
                        _xmin + dx_map,
                        _xmax + dx_map,
                        _ymin + dy_map,
                        _ymax + dy_map,
                    );

                    last_mouse_pos = Some((mx, my));
                }
            }
        } else {
            is_panning = false;
            last_mouse_pos = None;
        }

        // اسکرول برای زوم در حالت زوم
        if mode == "zoom" {
            if let Some(scroll) = window.get_scroll_wheel() {
                if let Some((mx, my)) = window.get_mouse_pos(MouseMode::Pass) {
                    let (xmin, xmax, ymin, ymax) = view_bbox;

                    let mx_rel = mx as f64 / WIDTH as f64;
                    let my_rel = my as f64 / HEIGHT as f64;

                    let zoom_factor = if scroll.1 > 0.0 { 0.9 } else { 1.1 };

                    let width = xmax - xmin;
                    let height = ymax - ymin;

                    let new_width = width * zoom_factor;
                    let new_height = height * zoom_factor;

                    let center_x = xmin + mx_rel * width;
                    let center_y = ymin + (1.0 - my_rel) * height;

                    view_bbox = (
                        center_x - new_width * mx_rel,
                        center_x + new_width * (1.0 - mx_rel),
                        center_y - new_height * (1.0 - my_rel),
                        center_y + new_height * my_rel,
                    );
                }
            }
        }

        // کلیک برای identify
        if mode == "identify" {
            if let Some((mx, my)) = window.get_mouse_pos(MouseMode::Pass) {
                if window.get_mouse_down(MouseButton::Left) {
                    let (_xmin, _xmax, _ymin, _ymax) = view_bbox;
                    let wx = _xmin + (mx as f64 / WIDTH as f64) * (_xmax - _xmin);
                    let wy = _ymin + (1.0 - (my as f64 / HEIGHT as f64)) * (_ymax - _ymin);

                    let mut closest_idx = None;
                    let mut closest_dist = f64::MAX;

                    for (i, shape) in shapes.iter().enumerate() {
                        let dist = match shape {
                            Shape::Point(x, y) => ((x - wx).powi(2) + (y - wy).powi(2)).sqrt(),
                            Shape::Poly(pts) => pts.iter()
                                .map(|(x, y)| ((x - wx).powi(2) + (y - wy).powi(2)).sqrt())
                                .fold(f64::MAX, |a, b| a.min(b)),
                        };
                        if dist < closest_dist && dist < 0.01 * (_xmax - _xmin) {
                            closest_dist = dist;
                            closest_idx = Some(i);
                        }
                    }

                    if let Some(idx) = closest_idx {
                        selected_shape_idx = Some(idx);
                        println!("📌 Selected index: {}", idx);
                        if let Some(rec) = dbf_records.get(idx) {
                            println!("🧾 DBF Record:");
                            for (k, v) in rec {
                                println!("  {}: {}", k, v);
                            }
                        }
                    } else {
                        selected_shape_idx = None;
                    }
                }
            }
        }

        buffer.fill(0xFFFFFF);
        let (_xmin, _xmax, _ymin, _ymax) = view_bbox;

        for (i, shape) in shapes.iter().enumerate() {
            let color = if Some(i) == selected_shape_idx { 0xFF0000 } else { 0x0000FF };
            match shape {
                Shape::Point(x, y) => {
                    let (px, py) = normalize(*x, *y, view_bbox);
                    draw_point(&mut buffer, px, py, color);
                }
                Shape::Poly(pts) => {
                    for w in pts.windows(2) {
                        let (x0, y0) = normalize(w[0].0, w[0].1, view_bbox);
                        let (x1, y1) = normalize(w[1].0, w[1].1, view_bbox);
                        draw_line(&mut buffer, x0, y0, x1, y1, color);
                    }
                }
            }
        }

        window.update_with_buffer(&buffer, WIDTH, HEIGHT)?;
    }

    Ok(())
}
