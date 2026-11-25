"""把 journal 版单图拼成一张成图（180 mm 宽，300 dpi）与一张补充图。

布局：
  行1  A fig1(56)  B fig2(56)  C fig3(56)
  行2  D fig4(85)  E fig5(85)
  行3  F fig6(85)  G fig9(85)
  行4  H fig8(85)  I fig11(85)
补充图  fig7(85) + fig10(85)
"""
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "figures" / "journal"
DPI = 300
MM = DPI / 25.4          # px per mm
PAGE_W_MM = 180.0
GAP_MM = 5.0
ROW_GAP_MM = 6.0
LETTER_MM = 5.0


def mm(v):
    return int(round(v * MM))


LAYOUT = [
    [("A", "fig1_pIC50_distribution.png", 56.0),
     ("B", "fig2_observed_vs_predicted.png", 56.0),
     ("C", "fig3_y_randomization.png", 56.0)],
    [("D", "fig4_william_plot.png", 85.0),
     ("E", "fig5_descriptor_importance.png", 85.0)],
    [("F", "fig6_query_predictions.png", 85.0),
     ("G", "fig9_le_vs_pIC50.png", 85.0)],
    [("H", "fig8_ligand_efficiency.png", 85.0),
     ("I", "fig11_dft_orbitals.png", 85.0)],
]
SUPP = [("S1", "fig7_docking_score.png", 85.0), ("S2", "fig10_structures.png", 85.0)]


def font(size_mm):
    px = int(size_mm * MM)
    for f in ("C:/Windows/Fonts/arialbd.ttf", "C:/Windows/Fonts/arial.ttf"):
        if Path(f).exists():
            return ImageFont.truetype(f, px)
    return ImageFont.load_default()


def build(layout, out_name, page_w):
    prepared = []
    for row in layout:
        items = []
        for letter, fname, w_mm in row:
            im = Image.open(SRC / fname).convert("RGB")
            w = mm(w_mm)
            h = int(round(im.height * w / im.width))
            items.append((letter, im.resize((w, h), Image.LANCZOS), w, h))
        nums = len(items)
        total_w = sum(i[2] for i in items) + mm(GAP_MM) * (nums - 1)
        prepared.append((items, total_w))

    canvas_w = mm(page_w)
    total_h = mm(LETTER_MM) * len(prepared) + sum(max(i[3] for i in r[0]) for r in prepared) \
        + mm(ROW_GAP_MM) * (len(prepared) - 1)
    canvas = Image.new("RGB", (canvas_w, total_h), "white")
    draw = ImageDraw.Draw(canvas)
    f = font(3.4)

    y = 0
    for items, total_w in prepared:
        x = (canvas_w - total_w) // 2
        row_h = max(i[3] for i in items)
        for letter, im, w, h in items:
            draw.text((x + 1, y), letter, fill="black", font=f)
            canvas.paste(im, (x, y + mm(LETTER_MM)))
            x += w + mm(GAP_MM)
        y += mm(LETTER_MM) + row_h + mm(ROW_GAP_MM)

    out = SRC / out_name
    canvas.save(out, dpi=(DPI, DPI))
    print(f"{out.name}: {canvas.width}x{canvas.height} px = "
          f"{canvas.width/MM:.0f}x{canvas.height/MM:.0f} mm")
    return out


build(LAYOUT, "Figure1-2_ODC1_QSAR_combined.png", PAGE_W_MM)
build([SUPP], "Figure_S1_docking_and_structures.png", PAGE_W_MM)
