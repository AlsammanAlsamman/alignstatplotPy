"""Regenerate the desktop-GUI screenshot used in the README.

Run from the repo root: `python packaging/make_screenshots.py`
(needs the `gui` extra installed, and a real Qt platform -- not
QT_QPA_PLATFORM=offscreen, which renders text as tofu boxes for the
native widgets here, even though the embedded Matplotlib canvas is fine).
"""
from PySide6.QtWidgets import QApplication

from alignstatplot_py.gui_qt.app import STYLE_SHEET, MainWindow
from alignstatplot_py.pipeline import run_pipeline

app = QApplication([])
app.setStyleSheet(STYLE_SHEET)
window = MainWindow()
window.resize(1280, 860)
window.fasta_field.setText("examples/Example_Small.fasta")
window.anno_field.setText("examples/Example_Small.anno")

result = run_pipeline(
    fasta_path="examples/Example_Small.fasta",
    anno_path="examples/Example_Small.anno",
    outdir="build/_screenshot_output",
    n_clusters=4,
    verbose=False,
    save_figures=False,
)
window._on_done(result)

for i in range(window.tabs.count()):
    if "circle" in window.tabs.tabText(i).lower():
        window.tabs.setCurrentIndex(i)
        break

window.show()
app.processEvents()
window.grab().save("assets/screenshot-desktop-gui.png")
print("wrote assets/screenshot-desktop-gui.png")
