from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
import argparse
import sys

def load_font(size: int = 16):
    """Try to load a truetype font, fallback to default."""
    try:
        return ImageFont.truetype("DejaVuSans.ttf", size)
    except IOError:
        return ImageFont.load_default()

def add_label_below_image(image: Image.Image, label: str, font: ImageFont.ImageFont) -> Image.Image:
    """Add a label below the image and return the new image."""
    # Get text size using getbbox (reliable for modern Pillow)
    bbox = font.getbbox(label)
    text_width = bbox[2] - bbox[0]
    text_height = bbox[3] - bbox[1]

    padding = 10
    new_img = Image.new("RGB", (image.width, image.height + text_height + padding), color="white")
    new_img.paste(image, (0, 0))

    draw = ImageDraw.Draw(new_img)
    text_x = (image.width - text_width) // 2
    text_y = image.height + (padding // 2)
    draw.text((text_x, text_y), label, fill="black", font=font)

    return new_img

def join_images(image_path_a: Path, image_path_b: Path, label_a: str, label_b: str, font: ImageFont.ImageFont) -> Image.Image:
    """Join two labeled images side-by-side with a vertical separator."""
    img_a = Image.open(image_path_a)
    img_b = Image.open(image_path_b)

    # # Resize to same height
    # if img_a.height != img_b.height:
    #     common_height = min(img_a.height, img_b.height)
    #     img_a = img_a.resize((int(img_a.width * common_height / img_a.height), common_height))
    #     img_b = img_b.resize((int(img_b.width * common_height / img_b.height), common_height))

    # Add labels below each image
    labeled_a = add_label_below_image(img_a, label_a, font)
    labeled_b = add_label_below_image(img_b, label_b, font)

    # Dimensions
    separator_width = 1  # 1 pixel line
    total_width = labeled_a.width + labeled_b.width + separator_width
    max_height = max(labeled_a.height, labeled_b.height)

    # Create new canvas with separator
    combined = Image.new("RGB", (total_width, max_height), color="white")
    combined.paste(labeled_a, (0, 0))
    combined.paste(labeled_b, (labeled_a.width + separator_width, 0))

    # Draw vertical separator line
    draw = ImageDraw.Draw(combined)
    line_x = labeled_a.width
    draw.line([(line_x, 0), (line_x, max_height)], fill="black", width=separator_width)

    return combined

def main(folder_a: Path, folder_b: Path, output_folder: Path):
    if not folder_a.is_dir() or not folder_b.is_dir():
        print("[ERROR] One or both input folders are invalid.", file=sys.stderr)
        return

    output_folder.mkdir(parents=True, exist_ok=True)
    font = load_font(40)

    files_a = {f.name for f in folder_a.iterdir() if f.is_file()}
    files_b = {f.name for f in folder_b.iterdir() if f.is_file()}
    common_files = files_a & files_b

    missing_in_a = files_b - files_a
    missing_in_b = files_a - files_b

    for missing in sorted(missing_in_a):
        print(f"[ERROR] File '{missing}' is in Folder B but missing in Folder A", file=sys.stderr)
    for missing in sorted(missing_in_b):
        print(f"[ERROR] File '{missing}' is in Folder A but missing in Folder B", file=sys.stderr)

    label_a = folder_a.name
    label_b = folder_b.name

    for filename in sorted(common_files):
        try:
            path_a = folder_a / filename
            path_b = folder_b / filename
            output_path = output_folder / filename

            new_img = join_images(path_a, path_b, label_a, label_b, font)
            new_img.save(output_path)
            print(f"[INFO] Joined and saved: {output_path}")
        except Exception as e:
            print(f"[ERROR] Failed to process '{filename}': {e}", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Join matching images from two folders side-by-side with labels.")
    parser.add_argument("folderA", type=Path, help="Path to Folder A")
    parser.add_argument("folderB", type=Path, help="Path to Folder B")
    parser.add_argument("output_folder", type=Path, help="Path to output folder")

    args = parser.parse_args()
    main(args.folderA, args.folderB, args.output_folder)
