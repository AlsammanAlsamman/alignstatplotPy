"""Desktop GUI (PySide6): pick a FASTA (+ optional annotation) file, set
parameters, run the same :func:`alignstatplot_py.pipeline.run_pipeline`
used by the CLI and web GUI, and browse the resulting plots/tables.
"""
from __future__ import annotations

import sys
import traceback
from pathlib import Path

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from PySide6.QtCore import Qt, QThread, Signal
from PySide6.QtGui import QIcon
from PySide6.QtWidgets import (
    QApplication,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from ..pipeline import run_pipeline

ICON_PATH = Path(__file__).resolve().parent.parent / "data" / "icon.png"

STYLE_SHEET = """
QMainWindow { background-color: #FAFAFA; }
QPushButton {
    background-color: #5E35B1; color: white; border-radius: 6px;
    padding: 8px 16px; font-weight: 600;
}
QPushButton:hover { background-color: #4527A0; }
QPushButton:disabled { background-color: #B0A8C7; }
QLabel#Header { font-size: 20px; font-weight: 700; color: #311B92; }
QTabWidget::pane { border: 1px solid #E0E0E0; }
"""


class PipelineWorker(QThread):
    finished_ok = Signal(object)
    failed = Signal(str)

    def __init__(self, kwargs):
        super().__init__()
        self.kwargs = kwargs

    def run(self):
        try:
            result = run_pipeline(**self.kwargs)
            self.finished_ok.emit(result)
        except Exception:
            self.failed.emit(traceback.format_exc())


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("alignstatplot")
        if ICON_PATH.exists():
            self.setWindowIcon(QIcon(str(ICON_PATH)))
        self.resize(1200, 800)
        self.fasta_path = ""
        self.anno_path = ""
        self.worker = None

        root = QWidget()
        self.setCentralWidget(root)
        layout = QVBoxLayout(root)

        header = QLabel("alignstatplot")
        header.setObjectName("Header")
        layout.addWidget(header)
        layout.addWidget(QLabel("Sequence alignment statistics and fancy plots"))

        form_row = QHBoxLayout()
        layout.addLayout(form_row)

        form = QFormLayout()
        form_row.addLayout(form, stretch=2)

        self.fasta_field = QLineEdit()
        self.fasta_field.setPlaceholderText("Choose a FASTA file...")
        fasta_btn = QPushButton("Browse")
        fasta_btn.clicked.connect(self._choose_fasta)
        fasta_row = QHBoxLayout()
        fasta_row.addWidget(self.fasta_field)
        fasta_row.addWidget(fasta_btn)
        form.addRow("FASTA file:", fasta_row)

        self.anno_field = QLineEdit()
        self.anno_field.setPlaceholderText("Optional annotation file...")
        anno_btn = QPushButton("Browse")
        anno_btn.clicked.connect(self._choose_anno)
        anno_row = QHBoxLayout()
        anno_row.addWidget(self.anno_field)
        anno_row.addWidget(anno_btn)
        form.addRow("Annotation file:", anno_row)

        self.outdir_field = QLineEdit("output")
        form.addRow("Output folder:", self.outdir_field)

        self.max_miss = QDoubleSpinBox()
        self.max_miss.setRange(0.0, 1.0)
        self.max_miss.setSingleStep(0.05)
        self.max_miss.setValue(0.2)
        form.addRow("Max missing per column:", self.max_miss)

        self.n_clusters = QSpinBox()
        self.n_clusters.setRange(1, 20)
        self.n_clusters.setValue(4)
        form.addRow("Number of SNP clusters:", self.n_clusters)

        self.min_cluster_len = QSpinBox()
        self.min_cluster_len.setRange(1, 50)
        self.min_cluster_len.setValue(3)
        form.addRow("Min SNP cluster length:", self.min_cluster_len)

        run_col = QVBoxLayout()
        form_row.addLayout(run_col, stretch=1)
        self.run_btn = QPushButton("Run analysis")
        self.run_btn.clicked.connect(self._run_pipeline)
        run_col.addWidget(self.run_btn)
        self.export_btn = QPushButton("Export all plots...")
        self.export_btn.clicked.connect(self._export_plots)
        self.export_btn.setEnabled(False)
        run_col.addWidget(self.export_btn)
        run_col.addStretch()

        self.log_box = QTextEdit()
        self.log_box.setReadOnly(True)
        self.log_box.setMaximumHeight(90)
        layout.addWidget(self.log_box)

        self.tabs = QTabWidget()
        layout.addWidget(self.tabs, stretch=1)

        self._result = None

    def _choose_fasta(self):
        path, _ = QFileDialog.getOpenFileName(self, "Choose FASTA file", "", "FASTA (*.fasta *.fa *.fna);;All files (*)")
        if path:
            self.fasta_field.setText(path)

    def _choose_anno(self):
        path, _ = QFileDialog.getOpenFileName(self, "Choose annotation file", "", "All files (*)")
        if path:
            self.anno_field.setText(path)

    def _log(self, msg: str):
        self.log_box.append(msg)

    def _run_pipeline(self):
        fasta = self.fasta_field.text().strip()
        if not fasta:
            QMessageBox.warning(self, "Missing input", "Please choose a FASTA file first.")
            return
        anno = self.anno_field.text().strip() or None
        outdir = self.outdir_field.text().strip() or "output"

        self.run_btn.setEnabled(False)
        self.log_box.clear()
        self._log("Running analysis...")

        kwargs = dict(
            fasta_path=fasta,
            anno_path=anno,
            outdir=outdir,
            max_miss_per=self.max_miss.value(),
            n_clusters=self.n_clusters.value(),
            min_cluster_length=self.min_cluster_len.value(),
            verbose=True,
        )
        self.worker = PipelineWorker(kwargs)
        self.worker.finished_ok.connect(self._on_done)
        self.worker.failed.connect(self._on_failed)
        self.worker.start()

    def _on_failed(self, message: str):
        self.run_btn.setEnabled(False)
        self.run_btn.setEnabled(True)
        self._log("Failed:\n" + message)
        QMessageBox.critical(self, "Analysis failed", message)

    def _on_done(self, result):
        self._result = result
        self.run_btn.setEnabled(True)
        self.export_btn.setEnabled(True)
        self._log("Done.")
        self._populate_tabs(result)

    def _populate_tabs(self, result):
        self.tabs.clear()

        stats_tab = self._table_widget(result.stats_df)
        self.tabs.addTab(stats_tab, "Sequence stats")

        div_tab = self._table_widget(result.diversity_df)
        self.tabs.addTab(div_tab, "Diversity")

        for name, fig in result.figures.items():
            scroll = QScrollArea()
            canvas = FigureCanvasQTAgg(fig)
            scroll.setWidget(canvas)
            scroll.setWidgetResizable(True)
            self.tabs.addTab(scroll, name.replace("_", " ").title())

    def _table_widget(self, df) -> QTableWidget:
        table = QTableWidget(len(df), len(df.columns))
        table.setHorizontalHeaderLabels([str(c) for c in df.columns])
        for i, (_, row) in enumerate(df.iterrows()):
            for j, val in enumerate(row):
                table.setItem(i, j, QTableWidgetItem(str(val)))
        return table

    def _export_plots(self):
        if self._result is None:
            return
        directory = QFileDialog.getExistingDirectory(self, "Choose export folder")
        if not directory:
            return
        for name, fig in self._result.figures.items():
            fig.savefig(Path(directory) / f"{name}.png", dpi=200, bbox_inches="tight")
        self._result.stats_df.to_csv(Path(directory) / "sequence_stats.csv", index=False)
        self._result.diversity_df.to_csv(Path(directory) / "diversity_stats.csv", index=False)
        QMessageBox.information(self, "Export complete", f"Saved results to {directory}")


def main() -> None:
    app = QApplication(sys.argv)
    app.setStyleSheet(STYLE_SHEET)
    if ICON_PATH.exists():
        app.setWindowIcon(QIcon(str(ICON_PATH)))
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
