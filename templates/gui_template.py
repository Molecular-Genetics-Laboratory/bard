# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'qtcreator_out.ui'
#
# Created by: PyQt5 UI code generator 5.9.2
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(654, 763)
        MainWindow.setMinimumSize(QtCore.QSize(654, 763))
        MainWindow.setMaximumSize(QtCore.QSize(654, 763))
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.tabWidget = QtWidgets.QTabWidget(self.centralwidget)
        self.tabWidget.setGeometry(QtCore.QRect(20, 20, 611, 651))
        self.tabWidget.setObjectName("tabWidget")
        self.tab = QtWidgets.QWidget()
        self.tab.setObjectName("tab")
        self.lbl_annotation = QtWidgets.QLabel(self.tab)
        self.lbl_annotation.setGeometry(QtCore.QRect(30, 19, 181, 21))
        self.lbl_annotation.setObjectName("lbl_annotation")
        self.lbl_gtf_featuretag = QtWidgets.QLabel(self.tab)
        self.lbl_gtf_featuretag.setGeometry(QtCore.QRect(30, 61, 161, 16))
        self.lbl_gtf_featuretag.setObjectName("lbl_gtf_featuretag")
        self.lbl_gtf_featuretype = QtWidgets.QLabel(self.tab)
        self.lbl_gtf_featuretype.setGeometry(QtCore.QRect(30, 103, 171, 16))
        self.lbl_gtf_featuretype.setObjectName("lbl_gtf_featuretype")
        self.lbl_covcutoff = QtWidgets.QLabel(self.tab)
        self.lbl_covcutoff.setGeometry(QtCore.QRect(30, 287, 141, 16))
        self.lbl_covcutoff.setObjectName("lbl_covcutoff")
        self.lbl_offset_dir = QtWidgets.QLabel(self.tab)
        self.lbl_offset_dir.setGeometry(QtCore.QRect(30, 241, 131, 16))
        self.lbl_offset_dir.setObjectName("lbl_offset_dir")
        self.lbl_initscan = QtWidgets.QLabel(self.tab)
        self.lbl_initscan.setGeometry(QtCore.QRect(30, 371, 171, 21))
        self.lbl_initscan.setObjectName("lbl_initscan")
        self.lbl_overlap = QtWidgets.QLabel(self.tab)
        self.lbl_overlap.setGeometry(QtCore.QRect(30, 330, 161, 21))
        self.lbl_overlap.setObjectName("lbl_overlap")
        self.lbl_listaction = QtWidgets.QLabel(self.tab)
        self.lbl_listaction.setGeometry(QtCore.QRect(30, 492, 121, 21))
        self.lbl_listaction.setObjectName("lbl_listaction")
        self.lbl_kmer = QtWidgets.QLabel(self.tab)
        self.lbl_kmer.setGeometry(QtCore.QRect(30, 412, 131, 21))
        self.lbl_kmer.setObjectName("lbl_kmer")
        self.lbl_genelist = QtWidgets.QLabel(self.tab)
        self.lbl_genelist.setGeometry(QtCore.QRect(30, 451, 151, 21))
        self.lbl_genelist.setObjectName("lbl_genelist")
        self.lbl_overlapallow = QtWidgets.QLabel(self.tab)
        self.lbl_overlapallow.setGeometry(QtCore.QRect(30, 532, 181, 21))
        self.lbl_overlapallow.setObjectName("lbl_overlapallow")
        self.lbl_operonmember = QtWidgets.QLabel(self.tab)
        self.lbl_operonmember.setGeometry(QtCore.QRect(30, 571, 181, 21))
        self.lbl_operonmember.setObjectName("lbl_operonmember")
        self.gtf_entrybox = QtWidgets.QLineEdit(self.tab)
        self.gtf_entrybox.setGeometry(QtCore.QRect(230, 20, 251, 23))
        self.gtf_entrybox.setCursor(QtGui.QCursor(QtCore.Qt.IBeamCursor))
        self.gtf_entrybox.setObjectName("gtf_entrybox")
        self.gtf_genetag_entrybox = QtWidgets.QLineEdit(self.tab)
        self.gtf_genetag_entrybox.setGeometry(QtCore.QRect(230, 60, 251, 23))
        self.gtf_genetag_entrybox.setObjectName("gtf_genetag_entrybox")
        self.gtf_featuretype_entrybox = QtWidgets.QLineEdit(self.tab)
        self.gtf_featuretype_entrybox.setGeometry(QtCore.QRect(230, 100, 251, 23))
        self.gtf_featuretype_entrybox.setObjectName("gtf_featuretype_entrybox")
        self.cds_entrybox = QtWidgets.QLineEdit(self.tab)
        self.cds_entrybox.setGeometry(QtCore.QRect(230, 144, 251, 23))
        self.cds_entrybox.setObjectName("cds_entrybox")
        self.offset_rb_5 = QtWidgets.QRadioButton(self.tab)
        self.offset_rb_5.setGeometry(QtCore.QRect(231, 240, 100, 21))
        self.offset_rb_5.setObjectName("offset_rb_5")
        self.BG1 = QtWidgets.QButtonGroup(MainWindow)
        self.BG1.setObjectName("BG1")
        self.BG1.addButton(self.offset_rb_5)
        self.offset_rb_3 = QtWidgets.QRadioButton(self.tab)
        self.offset_rb_3.setGeometry(QtCore.QRect(380, 240, 100, 21))
        self.offset_rb_3.setObjectName("offset_rb_3")
        self.BG1.addButton(self.offset_rb_3)
        self.overlap_ignore_rb_yes = QtWidgets.QRadioButton(self.tab)
        self.overlap_ignore_rb_yes.setGeometry(QtCore.QRect(231, 330, 100, 21))
        self.overlap_ignore_rb_yes.setObjectName("overlap_ignore_rb_yes")
        self.buttonGroup = QtWidgets.QButtonGroup(MainWindow)
        self.buttonGroup.setObjectName("buttonGroup")
        self.buttonGroup.addButton(self.overlap_ignore_rb_yes)
        self.overlap_ignore_rb_no = QtWidgets.QRadioButton(self.tab)
        self.overlap_ignore_rb_no.setGeometry(QtCore.QRect(380, 330, 100, 21))
        self.overlap_ignore_rb_no.setObjectName("overlap_ignore_rb_no")
        self.buttonGroup.addButton(self.overlap_ignore_rb_no)
        self.initscan_lower_entrybox = QtWidgets.QLineEdit(self.tab)
        self.initscan_lower_entrybox.setGeometry(QtCore.QRect(230, 370, 101, 23))
        self.initscan_lower_entrybox.setObjectName("initscan_lower_entrybox")
        self.label_35 = QtWidgets.QLabel(self.tab)
        self.label_35.setGeometry(QtCore.QRect(350, 370, 31, 21))
        self.label_35.setObjectName("label_35")
        self.initscan_upper_entrybox = QtWidgets.QLineEdit(self.tab)
        self.initscan_upper_entrybox.setGeometry(QtCore.QRect(380, 370, 101, 23))
        self.initscan_upper_entrybox.setObjectName("initscan_upper_entrybox")
        self.kmer_entrybox = QtWidgets.QLineEdit(self.tab)
        self.kmer_entrybox.setGeometry(QtCore.QRect(230, 410, 251, 23))
        self.kmer_entrybox.setObjectName("kmer_entrybox")
        self.genelist_entrybox = QtWidgets.QLineEdit(self.tab)
        self.genelist_entrybox.setGeometry(QtCore.QRect(230, 450, 251, 23))
        self.genelist_entrybox.setObjectName("genelist_entrybox")
        self.btn_gtf = QtWidgets.QPushButton(self.tab)
        self.btn_gtf.setGeometry(QtCore.QRect(510, 20, 80, 23))
        self.btn_gtf.setObjectName("btn_gtf")
        self.btn_cds = QtWidgets.QPushButton(self.tab)
        self.btn_cds.setGeometry(QtCore.QRect(510, 144, 80, 23))
        self.btn_cds.setObjectName("btn_cds")
        self.btn_genelist = QtWidgets.QPushButton(self.tab)
        self.btn_genelist.setGeometry(QtCore.QRect(510, 450, 80, 23))
        self.btn_genelist.setObjectName("btn_genelist")
        self.lbl_cds_file = QtWidgets.QLabel(self.tab)
        self.lbl_cds_file.setGeometry(QtCore.QRect(30, 148, 171, 16))
        self.lbl_cds_file.setObjectName("lbl_cds_file")
        self.cov_metric_dropdown = QtWidgets.QComboBox(self.tab)
        self.cov_metric_dropdown.setGeometry(QtCore.QRect(368, 282, 111, 23))
        self.cov_metric_dropdown.setObjectName("cov_metric_dropdown")
        self.cov_cutoff_spinbox = QtWidgets.QDoubleSpinBox(self.tab)
        self.cov_cutoff_spinbox.setGeometry(QtCore.QRect(230, 282, 111, 24))
        self.cov_cutoff_spinbox.setObjectName("cov_cutoff_spinbox")
        self.listaction_dropdown = QtWidgets.QComboBox(self.tab)
        self.listaction_dropdown.setGeometry(QtCore.QRect(230, 490, 251, 23))
        self.listaction_dropdown.setObjectName("listaction_dropdown")
        self.allowoverlap_entrybox = QtWidgets.QLineEdit(self.tab)
        self.allowoverlap_entrybox.setGeometry(QtCore.QRect(230, 530, 251, 23))
        self.allowoverlap_entrybox.setObjectName("allowoverlap_entrybox")
        self.btn_overlap = QtWidgets.QPushButton(self.tab)
        self.btn_overlap.setGeometry(QtCore.QRect(510, 530, 80, 23))
        self.btn_overlap.setObjectName("btn_overlap")
        self.operonmember_entrybox = QtWidgets.QLineEdit(self.tab)
        self.operonmember_entrybox.setGeometry(QtCore.QRect(230, 570, 251, 23))
        self.operonmember_entrybox.setObjectName("operonmember_entrybox")
        self.btn_operon = QtWidgets.QPushButton(self.tab)
        self.btn_operon.setGeometry(QtCore.QRect(510, 570, 80, 23))
        self.btn_operon.setObjectName("btn_operon")
        self.bam_entrybox = QtWidgets.QLineEdit(self.tab)
        self.bam_entrybox.setGeometry(QtCore.QRect(230, 190, 251, 23))
        self.bam_entrybox.setObjectName("bam_entrybox")
        self.lbl_bamfile = QtWidgets.QLabel(self.tab)
        self.lbl_bamfile.setGeometry(QtCore.QRect(30, 193, 171, 16))
        self.lbl_bamfile.setObjectName("lbl_bamfile")
        self.btn_bam = QtWidgets.QPushButton(self.tab)
        self.btn_bam.setGeometry(QtCore.QRect(510, 190, 80, 23))
        self.btn_bam.setObjectName("btn_bam")
        self.tabWidget.addTab(self.tab, "")
        self.tab_2 = QtWidgets.QWidget()
        self.tab_2.setObjectName("tab_2")
        self.label = QtWidgets.QLabel(self.tab_2)
        self.label.setGeometry(QtCore.QRect(19, 4, 501, 41))
        self.label.setObjectName("label")
        self.lbl_ldtran_maxcut = QtWidgets.QLabel(self.tab_2)
        self.lbl_ldtran_maxcut.setGeometry(QtCore.QRect(48, 138, 160, 21))
        self.lbl_ldtran_maxcut.setObjectName("lbl_ldtran_maxcut")
        self.ldtran_sum_cutoff_spinbox = QtWidgets.QDoubleSpinBox(self.tab_2)
        self.ldtran_sum_cutoff_spinbox.setGeometry(QtCore.QRect(258, 184, 121, 24))
        self.ldtran_sum_cutoff_spinbox.setObjectName("ldtran_sum_cutoff_spinbox")
        self.lbl_ldtran_sumcut = QtWidgets.QLabel(self.tab_2)
        self.lbl_ldtran_sumcut.setGeometry(QtCore.QRect(48, 186, 181, 16))
        self.lbl_ldtran_sumcut.setObjectName("lbl_ldtran_sumcut")
        self.ldtran_max_cutoff_spinbox = QtWidgets.QDoubleSpinBox(self.tab_2)
        self.ldtran_max_cutoff_spinbox.setGeometry(QtCore.QRect(258, 140, 121, 24))
        self.ldtran_max_cutoff_spinbox.setObjectName("ldtran_max_cutoff_spinbox")
        self.lbl_rdcov_term = QtWidgets.QLabel(self.tab_2)
        self.lbl_rdcov_term.setGeometry(QtCore.QRect(48, 232, 181, 16))
        self.lbl_rdcov_term.setObjectName("lbl_rdcov_term")
        self.lbl_rdcov_nt_start = QtWidgets.QLabel(self.tab_2)
        self.lbl_rdcov_nt_start.setGeometry(QtCore.QRect(48, 277, 181, 16))
        self.lbl_rdcov_nt_start.setObjectName("lbl_rdcov_nt_start")
        self.lbl_verbose_filename = QtWidgets.QLabel(self.tab_2)
        self.lbl_verbose_filename.setGeometry(QtCore.QRect(48, 323, 181, 16))
        self.lbl_verbose_filename.setObjectName("lbl_verbose_filename")
        self.verbose_filename_checkbox = QtWidgets.QCheckBox(self.tab_2)
        self.verbose_filename_checkbox.setGeometry(QtCore.QRect(260, 320, 85, 21))
        self.verbose_filename_checkbox.setObjectName("verbose_filename_checkbox")
        self.dropdown_readcov_frm_terminal = QtWidgets.QComboBox(self.tab_2)
        self.dropdown_readcov_frm_terminal.setGeometry(QtCore.QRect(258, 229, 121, 24))
        self.dropdown_readcov_frm_terminal.setObjectName("dropdown_readcov_frm_terminal")
        self.gtf_transcriptTag_entrybox = QtWidgets.QLineEdit(self.tab_2)
        self.gtf_transcriptTag_entrybox.setGeometry(QtCore.QRect(258, 93, 121, 24))
        self.gtf_transcriptTag_entrybox.setObjectName("gtf_transcriptTag_entrybox")
        self.lbl_transcript_id = QtWidgets.QLabel(self.tab_2)
        self.lbl_transcript_id.setGeometry(QtCore.QRect(48, 91, 191, 21))
        self.lbl_transcript_id.setObjectName("lbl_transcript_id")
        self.ignore_nts_start_spinbox = QtWidgets.QSpinBox(self.tab_2)
        self.ignore_nts_start_spinbox.setGeometry(QtCore.QRect(258, 274, 121, 24))
        self.ignore_nts_start_spinbox.setObjectName("ignore_nts_start_spinbox")
        self.tabWidget.addTab(self.tab_2, "")
        self.buttonBox = QtWidgets.QDialogButtonBox(self.centralwidget)
        self.buttonBox.setGeometry(QtCore.QRect(460, 680, 171, 41))
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.btn_export_json = QtWidgets.QPushButton(self.centralwidget)
        self.btn_export_json.setGeometry(QtCore.QRect(36, 688, 131, 23))
        self.btn_export_json.setObjectName("btn_export_json")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 654, 20))
        self.menubar.setObjectName("menubar")
        self.menuFile = QtWidgets.QMenu(self.menubar)
        self.menuFile.setObjectName("menuFile")
        self.menuHelp = QtWidgets.QMenu(self.menubar)
        self.menuHelp.setObjectName("menuHelp")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)
        self.actionSave_configuration = QtWidgets.QAction(MainWindow)
        self.actionSave_configuration.setObjectName("actionSave_configuration")
        self.actionExit = QtWidgets.QAction(MainWindow)
        self.actionExit.setObjectName("actionExit")
        self.actionAbout = QtWidgets.QAction(MainWindow)
        self.actionAbout.setObjectName("actionAbout")
        self.actionLoad_configuration = QtWidgets.QAction(MainWindow)
        self.actionLoad_configuration.setObjectName("actionLoad_configuration")
        self.actionExit_2 = QtWidgets.QAction(MainWindow)
        self.actionExit_2.setObjectName("actionExit_2")
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionLoad_configuration)
        self.menuFile.addAction(self.actionExit_2)
        self.menuHelp.addAction(self.actionAbout)
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuHelp.menuAction())

        self.retranslateUi(MainWindow)
        self.tabWidget.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "bard - Configuration panel"))
        self.lbl_annotation.setText(_translate("MainWindow", "Annotation file (GFF/GTF) *"))
        self.lbl_gtf_featuretag.setText(_translate("MainWindow", "Gene identifier tag *"))
        self.lbl_gtf_featuretype.setText(_translate("MainWindow", "Annotation feature type  *"))
        self.lbl_covcutoff.setText(_translate("MainWindow", "Read coverage cutoff"))
        self.lbl_offset_dir.setText(_translate("MainWindow", "P-site offset direction"))
        self.lbl_initscan.setText(_translate("MainWindow", "Initiation scan window (nt)"))
        self.lbl_overlap.setText(_translate("MainWindow", "Ignore overlapping genes"))
        self.lbl_listaction.setText(_translate("MainWindow", "Gene list action"))
        self.lbl_kmer.setText(_translate("MainWindow", "Use readlengths (nt)"))
        self.lbl_genelist.setText(_translate("MainWindow", "Arbitrary gene list (TXT)"))
        self.lbl_overlapallow.setText(_translate("MainWindow", "Genes overlap allowed (TXT)"))
        self.lbl_operonmember.setText(_translate("MainWindow", "Operon member genes (TXT)"))
        self.offset_rb_5.setText(_translate("MainWindow", "5\' end"))
        self.offset_rb_3.setText(_translate("MainWindow", "3\' end"))
        self.overlap_ignore_rb_yes.setText(_translate("MainWindow", "Yes"))
        self.overlap_ignore_rb_no.setText(_translate("MainWindow", "No"))
        self.label_35.setText(_translate("MainWindow", "to"))
        self.btn_gtf.setText(_translate("MainWindow", "Select file"))
        self.btn_cds.setText(_translate("MainWindow", "Select file"))
        self.btn_genelist.setText(_translate("MainWindow", "Select file"))
        self.lbl_cds_file.setText(_translate("MainWindow", "CDS file (FASTA) *"))
        self.btn_overlap.setText(_translate("MainWindow", "Select file"))
        self.btn_operon.setText(_translate("MainWindow", "Select file"))
        self.lbl_bamfile.setText(_translate("MainWindow", "Alignment file (BAM) *"))
        self.btn_bam.setText(_translate("MainWindow", "Select file"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab), _translate("MainWindow", "General"))
        self.label.setText(_translate("MainWindow", "Some internal defaults. Changing these may have unintended side effects."))
        self.lbl_ldtran_maxcut.setText(_translate("MainWindow", "LDTRAN_MAX_CUTOFF"))
        self.lbl_ldtran_sumcut.setText(_translate("MainWindow", "LDTRAN_SUM_CUTOFF"))
        self.lbl_rdcov_term.setText(_translate("MainWindow", "RDCOV_FRM_TERMINAL"))
        self.lbl_rdcov_nt_start.setText(_translate("MainWindow", "RDCOV_IGNORE_NTS_START"))
        self.lbl_verbose_filename.setText(_translate("MainWindow", "VERBOSE_FILENAMES"))
        self.verbose_filename_checkbox.setText(_translate("MainWindow", "Yes"))
        self.lbl_transcript_id.setText(_translate("MainWindow", "TRANSCRIPT_IDENTIFIER_TAG"))
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2), _translate("MainWindow", "Advanced"))
        self.btn_export_json.setText(_translate("MainWindow", "Save configuration"))
        self.menuFile.setTitle(_translate("MainWindow", "File"))
        self.menuHelp.setTitle(_translate("MainWindow", "Help"))
        self.actionSave_configuration.setText(_translate("MainWindow", "Export to JSON"))
        self.actionExit.setText(_translate("MainWindow", "Exit"))
        self.actionAbout.setText(_translate("MainWindow", "About"))
        self.actionLoad_configuration.setText(_translate("MainWindow", "Load configuration"))
        self.actionExit_2.setText(_translate("MainWindow", "Exit"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
