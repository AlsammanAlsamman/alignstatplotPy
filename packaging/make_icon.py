"""Generate the app icon: a small circular-alignment glyph (echoing the
package's signature 'Circular alignment overview' plot) on a purple
rounded-square background, exported as a multi-resolution .ico plus a
standalone PNG for the README.
"""
from PIL import Image, ImageDraw

ACCENT = (94, 53, 177, 255)
ACCENT_DARK = (49, 27, 146, 255)
NUC_COLORS = [
    (46, 125, 50, 255),  # A - green
    (21, 101, 192, 255),  # C - blue
    (198, 40, 40, 255),  # T - red
    (249, 168, 37, 255),  # G - orange
]

SIZE = 512


def make_base_image() -> Image.Image:
    img = Image.new("RGBA", (SIZE, SIZE), (0, 0, 0, 0))
    draw = ImageDraw.Draw(img)

    pad = 18
    radius = 110
    draw.rounded_rectangle(
        [pad, pad, SIZE - pad, SIZE - pad], radius=radius, fill=ACCENT, outline=ACCENT_DARK, width=6
    )

    center = SIZE / 2
    outer_r = SIZE * 0.34
    inner_r = SIZE * 0.15
    n = len(NUC_COLORS)
    gap = 6
    for i, color in enumerate(NUC_COLORS):
        start = 360 * i / n + gap
        end = 360 * (i + 1) / n - gap
        bbox = [center - outer_r, center - outer_r, center + outer_r, center + outer_r]
        draw.pieslice(bbox, start=start, end=end, fill=color)

    inner_bbox = [center - inner_r, center - inner_r, center + inner_r, center + inner_r]
    draw.ellipse(inner_bbox, fill=(255, 255, 255, 255))

    ring_r = SIZE * 0.40
    ring_bbox = [center - ring_r, center - ring_r, center + ring_r, center + ring_r]
    draw.ellipse(ring_bbox, outline=(255, 255, 255, 230), width=5)

    return img


if __name__ == "__main__":
    base = make_base_image()
    base.save("assets/icon.png")
    base.save("assets/icon.ico", sizes=[(16, 16), (24, 24), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)])
    print("wrote assets/icon.png and assets/icon.ico")
