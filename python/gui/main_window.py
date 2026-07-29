"""
PyQt5 graphical interface — Fluorescence spectrum analysis.

Architecture:
  - FluorescenceApp  : main window with tabs
      - "Raw data" tab           : raw 375 / 405 nm spectra
      - "Fluorescence fit" tab   : spectral decomposition
  - FitWorker        : computation thread (background fit)

Automatic behavior:
  - As soon as a file is selected -> spectra displayed + auto fit
  - As soon as a fluorophore is checked/unchecked or correction toggled -> fit re-run
  - As soon as the λ range changes -> spectra and fit refreshed
  - No "run fit" button
"""

import sys
import os
import traceback

import numpy as np

try:
    from PyQt5.QtWidgets import (
        QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
        QLabel, QLineEdit, QPushButton, QCheckBox, QRadioButton,
        QGroupBox, QFileDialog, QMessageBox, QSplitter, QScrollArea,
        QDoubleSpinBox, QProgressBar, QFrame, QTabWidget,
    )
    from PyQt5.QtCore import Qt, QThread, pyqtSignal, QObject, QTimer
    PYQT5_AVAILABLE = True
except ImportError:
    PYQT5_AVAILABLE = False

try:
    from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
    from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False


# ── Component color palette ──────────────────────────────────────────────────
COMP_COLORS = {
    'FAD':      '#8B2FC9',
    'NADH':     '#D4A017',
    'FMN':      '#1A9970',
    'Lipo':     '#2E86DE',
    'PpIX_636': '#C0392B',
    'PpIX_620': '#E67E22',
}
COMP_LABELS = {
    'FAD':      'FAD',
    'NADH':     'NADH',
    'FMN':      'protein-bound FMN',
    'Lipo':     'Lipopigments',
    'PpIX_636': 'PpIX 636 nm',
    'PpIX_620': 'PpIX 620 nm',
}

# ── QSS stylesheet ───────────────────────────────────────────────────────────
APP_STYLE = """
QMainWindow, QWidget#ControlPanel {
    background-color: #F5F6FA;
}
QGroupBox {
    font-size: 12px;
    font-weight: bold;
    color: #2C3E50;
    border: 1.5px solid #BDC3C7;
    border-radius: 6px;
    margin-top: 10px;
    padding-top: 6px;
}
QGroupBox::title {
    subcontrol-origin: margin;
    subcontrol-position: top left;
    left: 10px;
    padding: 0 4px;
    background-color: #F5F6FA;
}
QCheckBox {
    font-size: 11px;
    spacing: 6px;
    color: #2C3E50;
}
QCheckBox::indicator {
    width: 15px; height: 15px;
    border-radius: 3px;
    border: 1.5px solid #95A5A6;
    background: white;
}
QCheckBox::indicator:checked {
    background: #3498DB;
    border-color: #2980B9;
    image: url(none);
}
QRadioButton {
    font-size: 11px;
    spacing: 6px;
    color: #2C3E50;
}
QDoubleSpinBox {
    font-size: 11px;
    border: 1px solid #BDC3C7;
    border-radius: 4px;
    padding: 3px 6px;
    background: white;
    min-width: 70px;
}
QDoubleSpinBox:focus { border-color: #3498DB; }
QLineEdit {
    font-size: 11px;
    border: 1px solid #BDC3C7;
    border-radius: 4px;
    padding: 3px 6px;
    background: white;
}
QLineEdit:focus { border-color: #3498DB; }
QPushButton {
    font-size: 11px;
    border: 1px solid #BDC3C7;
    border-radius: 4px;
    padding: 4px 10px;
    background: #ECF0F1;
    color: #2C3E50;
}
QPushButton:hover  { background: #D5DBDB; }
QPushButton:pressed { background: #BDC3C7; }
QLabel#SectionTitle {
    font-size: 13px;
    font-weight: bold;
    color: #2C3E50;
}
QLabel#StatusLabel {
    font-size: 10px;
    color: #7F8C8D;
}
QProgressBar {
    border: 1px solid #BDC3C7;
    border-radius: 3px;
    text-align: center;
    height: 14px;
    font-size: 10px;
}
QProgressBar::chunk { background: #3498DB; border-radius: 3px; }
QTabWidget::pane {
    border: 1px solid #BDC3C7;
    border-radius: 0px 4px 4px 4px;
    background: white;
}
QTabBar::tab {
    font-size: 11px;
    padding: 7px 20px;
    background: #ECF0F1;
    border: 1px solid #BDC3C7;
    border-bottom: none;
    border-top-left-radius: 4px;
    border-top-right-radius: 4px;
    color: #2C3E50;
    min-width: 150px;
}
QTabBar::tab:selected {
    background: #3498DB;
    color: white;
    font-weight: bold;
    border-color: #2980B9;
}
QTabBar::tab:hover:!selected { background: #D5DBDB; }
"""

BTN_FILE_STYLE = """
QPushButton {
    background: #3498DB; color: white; font-weight: bold;
    border-radius: 5px; padding: 6px 14px; font-size: 12px;
}
QPushButton:hover  { background: #2980B9; }
QPushButton:pressed { background: #1F618D; }
"""


# ── Worker thread ─────────────────────────────────────────────────────────────
class FitWorker(QObject):
    finished = pyqtSignal(dict)
    error    = pyqtSignal(str)
    progress = pyqtSignal(str)

    def __init__(self, params):
        super().__init__()
        self.params = params

    def run(self):
        try:
            p = self.params
            sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            from functions.load_fluo import load_fluo
            from functions.corrected_fluo import corrected_fluo_ls_wl_1

            self.progress.emit("Loading FAD / NADH…")
            fad_385,  fad_405,  _ = load_fluo(p['fad_path'],  p['wl_min'], p['wl_max'])
            nadh_385, nadh_405, _ = load_fluo(p['nadh_path'], p['wl_min'], p['wl_max'])
            fluorophores_385 = np.vstack([fad_385,  nadh_385])
            fluorophores_405 = np.vstack([fad_405,  nadh_405])

            # Loading optional reference spectra (replace the gaussians)
            ref_spectra_385, ref_spectra_405 = {}, {}
            for name, path_ref in p.get('ref_spectra_paths', {}).items():
                if path_ref and os.path.isfile(path_ref):
                    try:
                        self.progress.emit(f"Loading spectrum {name}…")
                        s385, s405, lam_ref = load_fluo(path_ref, p['wl_min'], p['wl_max'])
                        # Pass (lambda, values) to allow fine-grained interpolation
                        ref_spectra_385[name] = (lam_ref, s385)
                        ref_spectra_405[name] = (lam_ref, s405)
                    except Exception as e:
                        self.progress.emit(f"Warning: spectrum {name}: {e}")

            self.progress.emit("Optimization in progress…")
            result = corrected_fluo_ls_wl_1(
                path=p['data_path'], file=p['filename'],
                name_spectralon=p['name_spectralon'],
                is_small_spectralon=p['is_small_spectralon'],
                min_wl_reflectance=p['wl_min'], max_wl_reflectance=p['wl_max'],
                min_wl_fluo=p['wl_min'],        max_wl_fluo=p['wl_max'],
                fluorophores_385=fluorophores_385,
                fluorophores_405=fluorophores_405,
                apply_optical_correction=p['apply_optical_correction'],
                active_fluorophores=p['active_fluorophores'],
                spectralon_theo_path=p.get('spectralon_theo_path'),
                scale_path=p.get('scale_path'),
                ref_spectra_385=ref_spectra_385 or None,
                ref_spectra_405=ref_spectra_405 or None,
            )
            result['_params'] = p
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(f"{type(e).__name__}: {str(e)}\n\n{traceback.format_exc()}")


# ── Main window ───────────────────────────────────────────────────────────────
class FluorescenceApp(QMainWindow):

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Fluorescence spectrum analysis")
        self.setGeometry(50, 50, 1400, 850)
        self.setStyleSheet(APP_STYLE)

        self._current_file = None
        self._current_path = None
        self._fit_thread   = None
        self._fit_worker   = None
        self._pending_fit  = False

        self._project_root         = self._find_project_root()
        self._fad_path             = self._auto_find('flavine_fluo.mat')
        self._nadh_path            = self._auto_find('NADH_fluo_4.mat')
        self._spectralon_theo_path = self._auto_find('Reflectance_values_array.txt')
        self._scale_path           = self._auto_find('scale_new_theory_intralipids.mat')

        self._update_timer = QTimer()
        self._update_timer.setSingleShot(True)
        self._update_timer.setInterval(400)
        self._update_timer.timeout.connect(self._refresh_displays)

        self._fit_timer = QTimer()
        self._fit_timer.setSingleShot(True)
        self._fit_timer.setInterval(1200)
        self._fit_timer.timeout.connect(self._trigger_fit)

        self._setup_ui()
        self._connect_signals()
        self.statusBar().showMessage("Select a *fluo.mat file to start")

    # ── UI construction ───────────────────────────────────────────────────────

    def _setup_ui(self):
        central = QWidget()
        central.setObjectName("ControlPanel")
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setSpacing(4)
        root.setContentsMargins(6, 6, 6, 4)

        root.addLayout(self._build_file_bar())

        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setFrameShadow(QFrame.Sunken)
        root.addWidget(sep)

        splitter = QSplitter(Qt.Horizontal)
        splitter.addWidget(self._build_left_panel())

        self._tabs = QTabWidget()
        self._tabs.addTab(self._build_raw_tab(), "  Raw data  ")
        self._tabs.addTab(self._build_reflectance_tab(), "  Reflectance  ")
        self._tabs.addTab(self._build_fit_tab(), "  Fluorescence fit  ")
        splitter.addWidget(self._tabs)

        splitter.setSizes([360, 1040])
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        root.addWidget(splitter, stretch=1)

        self._progress_bar = QProgressBar()
        self._progress_bar.setRange(0, 0)
        self._progress_bar.setVisible(False)
        self._progress_bar.setMaximumHeight(10)
        root.addWidget(self._progress_bar)

    def _build_file_bar(self):
        bar = QHBoxLayout()
        bar.setSpacing(8)

        lbl = QLabel("File:")
        lbl.setStyleSheet("font-weight:bold; font-size:12px;")
        bar.addWidget(lbl)

        self._file_label = QLineEdit("No file selected")
        self._file_label.setReadOnly(True)
        self._file_label.setStyleSheet("font-size:11px; color:#555;")
        bar.addWidget(self._file_label, stretch=1)

        self._btn_browse = QPushButton("  Browse…")
        self._btn_browse.setStyleSheet(BTN_FILE_STYLE)
        self._btn_browse.setMinimumWidth(120)
        bar.addWidget(self._btn_browse)

        return bar

    def _build_left_panel(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setMinimumWidth(310)
        scroll.setMaximumWidth(400)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setStyleSheet("QScrollArea{border:none; background:#F5F6FA;}")

        w = QWidget()
        w.setObjectName("ControlPanel")
        lay = QVBoxLayout(w)
        lay.setSpacing(10)
        lay.setContentsMargins(8, 8, 8, 8)

        lay.addWidget(self._build_display_group())
        lay.addWidget(self._build_wl_group())
        lay.addWidget(self._build_optical_group())
        lay.addWidget(self._build_fit_group())
        lay.addWidget(self._build_ref_files_group())
        lay.addStretch()

        status_frame = QFrame()
        status_frame.setFrameShape(QFrame.StyledPanel)
        status_frame.setStyleSheet(
            "QFrame{background:#EBF5FB;border-radius:5px;border:1px solid #AED6F1;}")
        sf_lay = QVBoxLayout(status_frame)
        sf_lay.setContentsMargins(8, 6, 8, 6)
        sf_lay.setSpacing(2)
        self._status_lbl = QLabel("—")
        self._status_lbl.setObjectName("StatusLabel")
        self._status_lbl.setWordWrap(True)
        self._status_lbl.setAlignment(Qt.AlignCenter)
        sf_lay.addWidget(self._status_lbl)
        lay.addWidget(status_frame)

        scroll.setWidget(w)
        return scroll

    def _build_display_group(self):
        g = QGroupBox("Spectrum display")
        lay = QVBoxLayout(g)
        lay.setSpacing(5)
        self._radio_mean = QRadioButton("Average of acquisitions")
        self._radio_all  = QRadioButton("Each individual acquisition")
        self._radio_mean.setChecked(True)
        lay.addWidget(self._radio_mean)
        lay.addWidget(self._radio_all)
        return g

    def _build_wl_group(self):
        g = QGroupBox("Wavelength range")
        lay = QHBoxLayout(g)
        lay.setSpacing(8)
        lay.addWidget(QLabel("Min:"))
        self._wl_min = QDoubleSpinBox()
        self._wl_min.setRange(300, 1000)
        self._wl_min.setValue(480)
        self._wl_min.setDecimals(0)
        self._wl_min.setSuffix(" nm")
        lay.addWidget(self._wl_min)
        lay.addSpacing(8)
        lay.addWidget(QLabel("Max:"))
        self._wl_max = QDoubleSpinBox()
        self._wl_max.setRange(300, 1000)
        self._wl_max.setValue(645)
        self._wl_max.setDecimals(0)
        self._wl_max.setSuffix(" nm")
        lay.addWidget(self._wl_max)
        return g

    def _build_optical_group(self):
        g = QGroupBox("Correction using optical properties")
        lay = QVBoxLayout(g)
        lay.setSpacing(5)
        self._cb_optical = QCheckBox("Enable correction (Kim model fit)")
        self._cb_optical.setToolTip(
            "Fits the Kim model to the diffuse reflectance.\n"
            "Requires an associated *diffuseReflectances.mat file.\n"
            "Without correction: R = 1 (original MATLAB behavior)."
        )
        lay.addWidget(self._cb_optical)
        note = QLabel("Without correction: R = 1  |  Raw data shown in the fit")
        note.setWordWrap(True)
        note.setStyleSheet("font-size:10px; color:#7F8C8D; margin-left:4px;")
        lay.addWidget(note)
        return g

    def _build_fit_group(self):
        g = QGroupBox("Fluorophores included in the fit")
        lay = QVBoxLayout(g)
        lay.setSpacing(4)
        self._cb_fluorophores = {}

        items = [
            ('FAD',      'FAD (Flavin Adenine Dinucleotide)'),
            ('NADH',     'NADH'),
            ('FMN',      'protein-bound FMN'),
            ('Lipo',     'Lipopigments'),
            ('PpIX_636', 'Protoporphyrin IX — 636 nm'),
            ('PpIX_620', 'Protoporphyrin IX — 620 nm'),
        ]
        for key, label in items:
            cb = QCheckBox(label)
            cb.setChecked(True)
            color = COMP_COLORS[key]
            cb.setStyleSheet(
                f"QCheckBox::indicator:checked{{background:{color};"
                f"border-color:{color};}}"
            )
            self._cb_fluorophores[key] = cb
            lay.addWidget(cb)

        return g

    def _build_ref_files_group(self):
        g = QGroupBox("Fluorophore reference files")
        lay = QVBoxLayout(g)
        lay.setSpacing(6)
        self._ref_spec_edits = {}

        def _ref_row(label, path, slot, placeholder="file not found", width=60):
            h = QHBoxLayout()
            h.setSpacing(4)
            lbl = QLabel(label)
            lbl.setFixedWidth(width)
            lbl.setStyleSheet("font-size:11px;")
            h.addWidget(lbl)
            edit = QLineEdit(path or "")
            edit.setPlaceholderText(placeholder)
            h.addWidget(edit, stretch=1)
            btn = QPushButton("…")
            btn.setFixedSize(26, 26)
            btn.clicked.connect(slot)
            h.addWidget(btn)
            return h, edit

        row_fad,  self._edit_fad  = _ref_row("FAD:",  self._fad_path,  self._browse_fad)
        row_nadh, self._edit_nadh = _ref_row("NADH:", self._nadh_path, self._browse_nadh)
        lay.addLayout(row_fad)
        lay.addLayout(row_nadh)

        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setFrameShadow(QFrame.Sunken)
        lay.addWidget(sep)

        note = QLabel("Optional spectra — replace the gaussians:")
        note.setStyleSheet("font-size:10px; color:#7F8C8D;")
        note.setWordWrap(True)
        lay.addWidget(note)

        gauss_items = [
            ('FMN',      'FMN:'),
            ('Lipo',     'Lipo:'),
            ('PpIX_636', 'PpIX 636:'),
            ('PpIX_620', 'PpIX 620:'),
        ]
        for key, label in gauss_items:
            h = QHBoxLayout()
            h.setSpacing(4)
            lbl = QLabel(label)
            lbl.setFixedWidth(60)
            lbl.setStyleSheet("font-size:11px;")
            h.addWidget(lbl)
            edit = QLineEdit()
            edit.setPlaceholderText("optional — otherwise gaussian")
            edit.setStyleSheet("font-size:10px;")
            h.addWidget(edit, stretch=1)
            btn = QPushButton("…")
            btn.setFixedSize(22, 22)
            btn.setToolTip(f"Choose a reference spectrum for {key}\n"
                           "(replaces the gaussian in the fit)")
            btn.clicked.connect(lambda _, k=key: self._browse_ref_spectrum(k))
            h.addWidget(btn)
            self._ref_spec_edits[key] = edit
            lay.addLayout(h)

        return g

    def _build_raw_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(0, 0, 0, 0)
        self._fig_raw, (self._ax385, self._ax405) = plt.subplots(
            1, 2, figsize=(14, 5), tight_layout=True
        )
        self._canvas_raw = FigureCanvas(self._fig_raw)
        self._tb_raw = NavigationToolbar(self._canvas_raw, w)
        lay.addWidget(self._tb_raw)
        lay.addWidget(self._canvas_raw)
        for ax, laser in [(self._ax385, "375 nm"), (self._ax405, "405 nm")]:
            ax.set_xlabel("Wavelength (nm)", fontsize=10)
            ax.set_ylabel("Intensity (a.u.)", fontsize=10)
            ax.set_title(f"Fluorescence — Excitation {laser}", fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3)
        return w

    def _build_reflectance_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.setSpacing(4)

        ctrl = QHBoxLayout()
        self._cb_calibrate_refl = QCheckBox("Calibrated reflectance (spectralon)")
        self._cb_calibrate_refl.setToolTip(
            "Divides the raw signal by the spectralon spectrum and the scale\n"
            "factor. Requires Reflectance_values_array.txt,\n"
            "scale_new_theory_intralipids.mat, and an experimental\n"
            "spectralon file in the data folder."
        )
        self._cb_small_spectralon = QCheckBox("Small spectralon (÷0.8)")
        ctrl.addWidget(self._cb_calibrate_refl)
        ctrl.addWidget(self._cb_small_spectralon)
        ctrl.addStretch()
        self._lbl_refl_status = QLabel("—")
        self._lbl_refl_status.setStyleSheet("font-size:10px; color:#7F8C8D;")
        self._lbl_refl_status.setWordWrap(True)
        ctrl.addWidget(self._lbl_refl_status, stretch=1)
        lay.addLayout(ctrl)

        self._fig_refl, (self._ax_refl_short, self._ax_refl_long) = plt.subplots(
            1, 2, figsize=(14, 5), tight_layout=True
        )
        self._canvas_refl = FigureCanvas(self._fig_refl)
        self._tb_refl = NavigationToolbar(self._canvas_refl, w)
        lay.addWidget(self._tb_refl)
        lay.addWidget(self._canvas_refl)
        for ax, label in [(self._ax_refl_short, "short wavelength (Short)"),
                          (self._ax_refl_long, "long wavelength (Long)")]:
            ax.set_xlabel("Wavelength (nm)", fontsize=10)
            ax.set_ylabel("Reflectance (a.u.)", fontsize=10)
            ax.set_title(f"Diffuse reflectance — {label}", fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3)
        return w

    def _build_fit_tab(self):
        w = QWidget()
        lay = QVBoxLayout(w)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.setSpacing(4)

        self._fig_fit = plt.figure(figsize=(18, 10))
        gs = gridspec.GridSpec(2, 2, figure=self._fig_fit,
                               height_ratios=[5, 1], hspace=0.35, wspace=0.25)
        self._ax385_fit  = self._fig_fit.add_subplot(gs[0, 0])
        self._ax405_fit  = self._fig_fit.add_subplot(gs[0, 1])
        self._ax375_hist = self._fig_fit.add_subplot(gs[1, 0])
        self._ax405_hist = self._fig_fit.add_subplot(gs[1, 1])

        self._canvas_fit = FigureCanvas(self._fig_fit)
        self._tb_fit = NavigationToolbar(self._canvas_fit, w)
        lay.addWidget(self._tb_fit)
        lay.addWidget(self._canvas_fit)

        # Info bar (alpha + redox)
        info_frame = QFrame()
        info_frame.setFrameShape(QFrame.StyledPanel)
        info_frame.setStyleSheet("background:#EBF5FB; border-radius:5px; padding:2px;")
        info_lay = QHBoxLayout(info_frame)
        info_lay.setContentsMargins(12, 4, 12, 4)

        self._lbl_alpha_fit = QLabel("α: —")
        self._lbl_alpha_fit.setStyleSheet(
            "font-size:12px; font-weight:bold; color:#1A5276;")
        self._lbl_redox_fit = QLabel("Redox: —")
        self._lbl_redox_fit.setStyleSheet(
            "font-size:12px; font-weight:bold; color:#1E8449;")
        self._lbl_fit_status = QLabel("")
        self._lbl_fit_status.setStyleSheet("font-size:11px; color:#7F8C8D;")
        self._lbl_fit_status.setAlignment(Qt.AlignRight)

        info_lay.addWidget(self._lbl_alpha_fit)
        info_lay.addSpacing(30)
        info_lay.addWidget(self._lbl_redox_fit)
        info_lay.addStretch()
        info_lay.addWidget(self._lbl_fit_status)
        lay.addWidget(info_frame)

        for ax, laser in [(self._ax385_fit, "375 nm"), (self._ax405_fit, "405 nm")]:
            ax.set_xlabel("Wavelength (nm)", fontsize=10)
            ax.set_ylabel("Normalized intensity (a.u.)", fontsize=10)
            ax.set_title(f"Fit — Excitation {laser}", fontsize=11, fontweight='bold')
            ax.grid(True, alpha=0.3)
        for ax_h in [self._ax375_hist, self._ax405_hist]:
            ax_h.set_ylabel("Fraction (%)", fontsize=9)
            ax_h.set_ylim(0, 110)
            ax_h.grid(True, alpha=0.3, axis='y')
        self._canvas_fit.draw()

        return w

    # ── Signal connections ────────────────────────────────────────────────────

    def _connect_signals(self):
        self._btn_browse.clicked.connect(self._browse_file)

        self._radio_mean.toggled.connect(self._schedule_display)
        self._radio_all.toggled.connect(self._schedule_display)
        self._wl_min.valueChanged.connect(self._schedule_display)
        self._wl_max.valueChanged.connect(self._schedule_display)
        self._cb_calibrate_refl.stateChanged.connect(self._schedule_display)
        self._cb_small_spectralon.stateChanged.connect(self._schedule_display)

        self._cb_optical.stateChanged.connect(self._schedule_fit)
        self._wl_min.valueChanged.connect(self._schedule_fit)
        self._wl_max.valueChanged.connect(self._schedule_fit)
        for cb in self._cb_fluorophores.values():
            cb.stateChanged.connect(self._schedule_fit)
        for edit in self._ref_spec_edits.values():
            edit.textChanged.connect(self._schedule_fit)

    def _schedule_display(self):
        if self._current_file:
            self._update_timer.start()

    def _refresh_displays(self):
        self._display_spectra()
        self._display_reflectance()

    def _schedule_fit(self):
        if self._current_file:
            self._fit_timer.start()

    # ── Files ─────────────────────────────────────────────────────────────────

    def _browse_file(self):
        start = self._current_path or self._project_root or ""
        path, _ = QFileDialog.getOpenFileName(
            self, "Select a fluorescence file",
            start, "Fluorescence files (*fluo.mat);;MATLAB files (*.mat);;All (*)"
        )
        if not path:
            return
        self._current_file = os.path.basename(path)
        self._current_path = os.path.dirname(path)
        self._file_label.setText(path)
        self.statusBar().showMessage(f"File: {self._current_file}")
        self._radio_mean.setChecked(True)
        self._display_spectra()
        self._display_reflectance()
        self._trigger_fit()

    def _browse_fad(self):
        p, _ = QFileDialog.getOpenFileName(
            self, "FAD — reference", self._current_path or "", "*.mat")
        if p:
            self._fad_path = p
            self._edit_fad.setText(p)

    def _browse_nadh(self):
        p, _ = QFileDialog.getOpenFileName(
            self, "NADH — reference", self._current_path or "", "*.mat")
        if p:
            self._nadh_path = p
            self._edit_nadh.setText(p)

    def _browse_ref_spectrum(self, fluorophore_name):
        p, _ = QFileDialog.getOpenFileName(
            self, f"Reference spectrum — {fluorophore_name}",
            self._current_path or "",
            "MATLAB files (*.mat);;Text files (*.txt);;All (*)"
        )
        if p:
            self._ref_spec_edits[fluorophore_name].setText(p)

    # ── Parameters ────────────────────────────────────────────────────────────

    def _get_params(self):
        return {
            'wl_min': float(self._wl_min.value()),
            'wl_max': float(self._wl_max.value()),
            'show_mean': self._radio_mean.isChecked(),
            'is_small_spectralon': False,
            'apply_optical_correction': self._cb_optical.isChecked(),
            'active_fluorophores': {k: cb.isChecked()
                                    for k, cb in self._cb_fluorophores.items()},
            'name_spectralon': 'spectralon_000_diffuseReflectances',
            'fad_path':  self._edit_fad.text().strip()  or self._fad_path  or "",
            'nadh_path': self._edit_nadh.text().strip() or self._nadh_path or "",
            'data_path': self._current_path,
            'filename':  self._current_file,
            'spectralon_theo_path': self._spectralon_theo_path,
            'scale_path': self._scale_path,
            'ref_spectra_paths': {k: e.text().strip()
                                  for k, e in self._ref_spec_edits.items()},
        }

    # ── Raw spectra display ───────────────────────────────────────────────────

    def _display_spectra(self):
        if not self._current_file:
            return
        p  = self._get_params()
        fp = os.path.join(self._current_path, self._current_file)
        try:
            sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            from functions.load_fluo import load_fluo
            all385, all405, lam = load_fluo(fp, p['wl_min'], p['wl_max'], return_all=True)
            self._plot_raw(lam, all385, all405, p['show_mean'])
            n = all385.shape[0]
            self.statusBar().showMessage(
                f"{self._current_file}  —  {n} acquisition(s)  —  "
                f"{p['wl_min']:.0f}–{p['wl_max']:.0f} nm"
            )
        except Exception as e:
            self.statusBar().showMessage(f"Error: {e}")

    def _plot_channel(self, ax, lam, data, color, title, ylabel, show_mean):
        """Plot the mean ± σ or each individual acquisition on an axis."""
        ax.clear()
        data = np.atleast_2d(data)
        n = data.shape[0]
        if show_mean:
            mu = np.mean(data, axis=0)
            ax.plot(lam, mu, color=color, lw=2, label='Mean')
            if n > 1:
                sig = np.std(data, axis=0)
                ax.fill_between(lam, mu - sig, mu + sig,
                                alpha=0.15, color=color, label='± 1 σ')
        else:
            cmap = plt.cm.viridis(np.linspace(0.1, 0.9, max(n, 1)))
            for i in range(n):
                ax.plot(lam, data[i], color=cmap[i], lw=1.2,
                        alpha=0.75, label=f'Acq. {i+1}')
        ax.set_xlabel("Wavelength (nm)", fontsize=10)
        ax.set_ylabel(ylabel, fontsize=10)
        ax.set_title(title, fontsize=11, fontweight='bold')
        ax.set_xlim(lam[0], lam[-1])
        ax.legend(loc='upper right', fontsize=8, framealpha=0.85,
                  edgecolor='#BDC3C7')
        ax.grid(True, alpha=0.3)

    def _plot_raw(self, lam, all385, all405, show_mean):
        self._plot_channel(self._ax385, lam, all385, '#2E86C1',
                           "Fluorescence — Excitation 375 nm", "Intensity (a.u.)", show_mean)
        self._plot_channel(self._ax405, lam, all405, '#C0392B',
                           "Fluorescence — Excitation 405 nm", "Intensity (a.u.)", show_mean)
        self._fig_raw.tight_layout()
        self._canvas_raw.draw()

    # ── Reflectance display (visualization only, no fit) ─────────────────────

    def _display_reflectance(self):
        if not self._current_file:
            return
        fp = os.path.join(self._current_path, self._current_file)
        wl_min = float(self._wl_min.value())
        wl_max = float(self._wl_max.value())
        show_mean = self._radio_mean.isChecked()
        calibrate = self._cb_calibrate_refl.isChecked()
        is_small = self._cb_small_spectralon.isChecked()

        try:
            sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
            from functions.reflectance_view import load_reflectance_view
            view = load_reflectance_view(
                fp,
                spectralon_theo_path=self._spectralon_theo_path,
                scale_path=self._scale_path,
                is_small_spectralon=is_small,
                wl_min=wl_min, wl_max=wl_max,
                return_all=True,
            )
        except FileNotFoundError as e:
            self._clear_reflectance_axes(str(e))
            return
        except Exception as e:
            self._clear_reflectance_axes(f"Loading error: {e}")
            return

        lam = view['lambda']
        if calibrate and view['R_Short'] is not None:
            data_short, data_long = view['R_Short'], view['R_Long']
            ylabel = "Calibrated reflectance (a.u.)"
        else:
            data_short, data_long = view['WLShort'], view['WLLong']
            ylabel = "Raw signal (a.u.)"

        self._plot_channel(self._ax_refl_short, lam, data_short, '#2E86C1',
                           "Reflectance — Short", ylabel, show_mean)
        self._plot_channel(self._ax_refl_long, lam, data_long, '#C0392B',
                           "Reflectance — Long", ylabel, show_mean)
        self._fig_refl.tight_layout()
        self._canvas_refl.draw()

        status = os.path.basename(view['reflectance_path'])
        if calibrate:
            if view['calibration_error']:
                status += f"  —  Calibration unavailable: {view['calibration_error']}  (raw signal shown)"
            else:
                status += "  —  calibrated reflectance"
        self._lbl_refl_status.setText(status)

    def _clear_reflectance_axes(self, message):
        for ax in (self._ax_refl_short, self._ax_refl_long):
            ax.clear()
            ax.text(0.5, 0.5, message, ha='center', va='center',
                    transform=ax.transAxes, fontsize=10, color='#7F8C8D', wrap=True)
            ax.set_xticks([])
            ax.set_yticks([])
        self._canvas_refl.draw()
        self._lbl_refl_status.setText(message)

    # ── Triggering the fit ────────────────────────────────────────────────────

    def _trigger_fit(self):
        if not self._current_file:
            return
        p = self._get_params()
        if not p['fad_path'] or not os.path.isfile(p['fad_path']):
            self._status_lbl.setText("FAD file not found — fit not possible.")
            return
        if not p['nadh_path'] or not os.path.isfile(p['nadh_path']):
            self._status_lbl.setText("NADH file not found — fit not possible.")
            return
        if p['apply_optical_correction']:
            refl = self._current_file.replace('fluo', 'diffuseReflectances')
            if not os.path.isfile(os.path.join(self._current_path, refl)):
                self._status_lbl.setText(
                    f"Reflectance file not found:\n{refl}\n"
                    "Disable the optical correction."
                )
                return

        if self._fit_thread and self._fit_thread.isRunning():
            self._pending_fit = True
            self._status_lbl.setText("Fit in progress… restart pending.")
            return

        self._start_fit(p)

    def _start_fit(self, p):
        self._pending_fit = False
        self._lbl_fit_status.setText("Optimization in progress…")
        self._progress_bar.setVisible(True)
        self._status_lbl.setText("Optimization in progress…")

        self._fit_thread = QThread()
        self._fit_worker = FitWorker(p)
        self._fit_worker.moveToThread(self._fit_thread)
        self._fit_thread.started.connect(self._fit_worker.run)
        self._fit_worker.finished.connect(self._on_fit_done)
        self._fit_worker.error.connect(self._on_fit_error)
        self._fit_worker.progress.connect(self._on_fit_progress)
        self._fit_worker.finished.connect(self._fit_thread.quit)
        self._fit_worker.error.connect(self._fit_thread.quit)
        self._fit_thread.finished.connect(self._on_thread_finished)
        self._fit_thread.start()

    def _on_fit_progress(self, msg):
        self._status_lbl.setText(msg)
        self._lbl_fit_status.setText(msg)
        self.statusBar().showMessage(msg)

    def _on_fit_done(self, result):
        self._update_fit_tab(result)
        self._lbl_fit_status.setText("")
        fluo = result.get('fluorophore', {})
        r385 = fluo.get('redox_385', 0.0)
        r405 = fluo.get('redox_405', 0.0)
        self._status_lbl.setText(
            f"Fit OK\nRedox 375 nm = {r385:.3f}\nRedox 405 nm = {r405:.3f}"
        )
        self.statusBar().showMessage(
            f"Fit complete  —  Redox 375 nm: {r385:.3f}  |  Redox 405 nm: {r405:.3f}"
        )
        self._tabs.setCurrentIndex(1)

    def _on_fit_error(self, msg):
        self._lbl_fit_status.setText("Fit error.")
        self._status_lbl.setText("Fit error.")
        self.statusBar().showMessage("Fit failed.")
        QMessageBox.critical(self, "Fit error", msg[:1500])

    def _on_thread_finished(self):
        self._progress_bar.setVisible(False)
        if self._pending_fit:
            self._trigger_fit()

    def _update_fit_tab(self, result: dict):
        p          = result.get('_params', {})
        apply_corr = p.get('apply_optical_correction', False)
        lam        = result['lambda_fluo']
        fluo       = result.get('fluorophore', {})
        res_385    = result.get('res_385')
        res_405    = result.get('res_405')
        alpha_385  = float(res_385[0]) if res_385 is not None else None
        alpha_405  = float(res_405[0]) if res_405 is not None else None

        name_map = {
            'NADH': 'NADH', 'FAD': 'flavine', 'FMN': 'gaussian',
            'Lipo': 'lipo', 'PpIX_620': 'PpIX_620', 'PpIX_636': 'PpIX_636',
        }
        hist_labels = {
            'FAD': 'FAD', 'NADH': 'NADH', 'FMN': 'FMN',
            'Lipo': 'Lipo', 'PpIX_636': 'PpIX\n636', 'PpIX_620': 'PpIX\n620',
        }

        for ax, ax_hist, suffix, raw_key, corr_key, laser in [
            (self._ax385_fit, self._ax375_hist,
             '_385_exp', 'S385total', 'S385_corrected', "375 nm"),
            (self._ax405_fit, self._ax405_hist,
             '_405_exp', 'S405total', 'S405_corrected', "405 nm"),
        ]:
            ax.clear()
            S_raw  = result.get(raw_key,  np.zeros_like(lam))
            S_corr = result.get(corr_key, S_raw)
            norm   = max(np.sum(S_corr), 1e-12)

            if not apply_corr:
                ax.plot(lam, S_raw / norm, color='#808B96', lw=1.5, ls='--',
                        alpha=0.7, label='Raw data (normalized)')
            else:
                ax.plot(lam, S_corr / norm, color='#2C3E50', lw=2,
                        label='Corrected signal', zorder=5)

            sum_fit = sum(
                fluo.get(name_map[k] + suffix, np.zeros_like(lam))
                for k in name_map
            )
            ax.plot(lam, sum_fit / norm, color='#2471A3', lw=2.5, alpha=0.85,
                    label='Sum of fit', zorder=4)

            active_comps = {}
            for k, mat_name in name_map.items():
                comp = fluo.get(mat_name + suffix, None)
                if comp is not None and np.any(comp != 0):
                    ax.plot(lam, comp / norm,
                            color=COMP_COLORS[k], lw=2, ls='--',
                            label=COMP_LABELS[k], zorder=3)
                    active_comps[k] = comp

            ax.set_xlabel("Wavelength (nm)", fontsize=10)
            ax.set_ylabel("Normalized intensity (a.u.)", fontsize=10)
            ax.set_title(f"Fit — Excitation {laser}", fontsize=11, fontweight='bold')
            ax.set_xlim(lam[0], lam[-1])
            ax.grid(True, alpha=0.3)
            ax.legend(loc='upper right', fontsize=8, framealpha=0.9,
                      edgecolor='#BDC3C7')

            # Histogram of relative fractions (integral of each component)
            ax_hist.clear()
            if active_comps:
                integrals = {k: max(float(np.trapz(v, lam)), 0.0)
                             for k, v in active_comps.items()}
                total = sum(integrals.values())
                if total > 0:
                    keys   = list(integrals.keys())
                    pcts   = [integrals[k] / total * 100 for k in keys]
                    colors = [COMP_COLORS[k] for k in keys]
                    xlbls  = [hist_labels[k] for k in keys]
                    bars = ax_hist.bar(range(len(keys)), pcts, color=colors, alpha=0.85,
                                      edgecolor='white', linewidth=0.5)
                    for bar, pct in zip(bars, pcts):
                        if pct >= 0.5:
                            ax_hist.text(bar.get_x() + bar.get_width() / 2,
                                         bar.get_height() + 1.0,
                                         f"{pct:.1f}%", ha='center', va='bottom',
                                         fontsize=7.5, color='#2C3E50')
                    ax_hist.set_xticks(range(len(keys)))
                    ax_hist.set_xticklabels(xlbls, fontsize=8)
            ax_hist.set_ylabel("Fraction (%)", fontsize=9)
            ax_hist.set_ylim(0, 115)
            ax_hist.grid(True, alpha=0.3, axis='y')
            ax_hist.set_title(f"Relative fractions — {laser}", fontsize=9)

        self._fig_fit.subplots_adjust(
            left=0.09, right=0.99, top=0.95, bottom=0.08,
            hspace=0.35, wspace=0.25
        )
        self._canvas_fit.draw()

        r385 = fluo.get('redox_385', None)
        r405 = fluo.get('redox_405', None)

        if apply_corr and alpha_385 is not None:
            self._lbl_alpha_fit.setText(
                f"α  —  375 nm: {alpha_385:.4f}   |   405 nm: {alpha_405:.4f}")
        else:
            self._lbl_alpha_fit.setText("No optical correction (α not applied)")

        if r385 is not None:
            self._lbl_redox_fit.setText(
                f"Redox  —  375 nm: {r385:.3f}   |   405 nm: {r405:.3f}")

    # ── Utilities ─────────────────────────────────────────────────────────────

    def _find_project_root(self):
        cur = os.path.dirname(os.path.abspath(__file__))
        for _ in range(6):
            if os.path.isdir(os.path.join(cur, 'matlab')) and \
               os.path.isdir(os.path.join(cur, 'data')):
                return cur
            parent = os.path.dirname(cur)
            if parent == cur:
                break
            cur = parent
        return None

    def _auto_find(self, filename):
        if not self._project_root:
            return None
        for root, dirs, files in os.walk(self._project_root):
            dirs[:] = [d for d in dirs if d not in ('.git', '__pycache__')]
            if filename in files:
                return os.path.join(root, filename)
        return None

    def closeEvent(self, event):
        if self._fit_thread and self._fit_thread.isRunning():
            self._fit_thread.quit()
            self._fit_thread.wait(3000)
        plt.close('all')
        event.accept()


# ── Entry point ───────────────────────────────────────────────────────────────

def run_fluorescence_app(argv=None):
    """
    Launches the application.  Importable from another script:
        from python.gui.main_window import run_fluorescence_app
        run_fluorescence_app()
    """
    if not PYQT5_AVAILABLE:
        raise ImportError("PyQt5 required: pip install PyQt5")
    if not MATPLOTLIB_AVAILABLE:
        raise ImportError("matplotlib required: pip install matplotlib")

    if argv is None:
        argv = sys.argv
    app = QApplication.instance()
    new_app = app is None
    if new_app:
        app = QApplication(argv)
        app.setStyle('Fusion')
    window = FluorescenceApp()
    window.show()
    if new_app:
        return app.exec_()
    return 0


if __name__ == '__main__':
    sys.exit(run_fluorescence_app())
