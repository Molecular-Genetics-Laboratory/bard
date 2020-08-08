from PyQt5 import QtCore, QtGui, QtWidgets
import json
import os

# Global configuration to be updated
global_config = {}

def save_json(data, file="object.json"):
    """
    Save object as JSON
    """
    with open(file, "w") as write_file:
        json.dump(data, write_file, indent=4)
        

def load_json(file="object.json"):
    """
    Load a JSON file
    """
    with open(file, "r") as read_file:
        return json.load(read_file)


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
        # -- [Set the default GTF attribute for a gene: ID/protein_id/etc
        self.gtf_genetag_entrybox.setText("ID")
        # -- ]
        self.gtf_featuretype_entrybox = QtWidgets.QLineEdit(self.tab)
        self.gtf_featuretype_entrybox.setGeometry(QtCore.QRect(230, 100, 251, 23))
        self.gtf_featuretype_entrybox.setObjectName("gtf_featuretype_entrybox")
        # -- [Set what type of feature we're looking at: gene/CDS/Exon etc
        self.gtf_featuretype_entrybox.setText("CDS")
        # -- ]
        self.cds_entrybox = QtWidgets.QLineEdit(self.tab)
        self.cds_entrybox.setGeometry(QtCore.QRect(230, 144, 251, 23))
        self.cds_entrybox.setObjectName("cds_entrybox")
        self.offset_rb_5 = QtWidgets.QRadioButton(self.tab)
        self.offset_rb_5.setGeometry(QtCore.QRect(231, 240, 100, 21))
        self.offset_rb_5.setObjectName("offset_rb_5")
        # -- [ Set 5' offset as default
        self.offset_rb_5.setChecked(True)
        # -- ]
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
        # -- [ Set ignore overlaps as the default
        self.overlap_ignore_rb_yes.setChecked(True)
        # -- ]
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
        # -- [ Range for scanning for initiation peak
        self.initscan_lower_entrybox.setText("-25")
        # --]
        self.label_35 = QtWidgets.QLabel(self.tab)
        self.label_35.setGeometry(QtCore.QRect(350, 370, 31, 21))
        self.label_35.setObjectName("label_35")
        self.initscan_upper_entrybox = QtWidgets.QLineEdit(self.tab)
        self.initscan_upper_entrybox.setGeometry(QtCore.QRect(380, 370, 101, 23))
        self.initscan_upper_entrybox.setObjectName("initscan_upper_entrybox")
        # -- [ Range for initiation scan window
        self.initscan_upper_entrybox.setText("-5")
        # -- ]
        self.kmer_entrybox = QtWidgets.QLineEdit(self.tab)
        self.kmer_entrybox.setGeometry(QtCore.QRect(230, 410, 251, 23))
        self.kmer_entrybox.setObjectName("kmer_entrybox")
        # -- [ Set readlengths to be used for the analysis
        self.kmer_entrybox.setText("25,26,27,28,29,30,31,32")
        # -- ]
        self.genelist_entrybox = QtWidgets.QLineEdit(self.tab)
        self.genelist_entrybox.setGeometry(QtCore.QRect(230, 450, 251, 23))
        self.genelist_entrybox.setObjectName("genelist_entrybox")
        self.btn_gtf = QtWidgets.QPushButton(self.tab)
        self.btn_gtf.setGeometry(QtCore.QRect(510, 20, 80, 23))
        self.btn_gtf.setObjectName("btn_gtf")
        # -- [Select annotation file
        self.btn_gtf.clicked.connect(lambda: self.choose_path_for_fileopen(self.gtf_entrybox))
        # -- ]
        self.btn_cds = QtWidgets.QPushButton(self.tab)
        self.btn_cds.setGeometry(QtCore.QRect(510, 144, 80, 23))
        self.btn_cds.setObjectName("btn_cds")
        # -- [Select CDS file
        self.btn_cds.clicked.connect(lambda: self.choose_path_for_fileopen(self.cds_entrybox))
        # -- ]
        self.btn_genelist = QtWidgets.QPushButton(self.tab)
        self.btn_genelist.setGeometry(QtCore.QRect(510, 450, 80, 23))
        self.btn_genelist.setObjectName("btn_genelist")
        # -- [ Select file for arbitrary gene list
        self.btn_genelist.clicked.connect(lambda: self.choose_path_for_fileopen(self.genelist_entrybox))
        # -- ]
        self.lbl_cds_file = QtWidgets.QLabel(self.tab)
        self.lbl_cds_file.setGeometry(QtCore.QRect(30, 148, 171, 16))
        self.lbl_cds_file.setObjectName("lbl_cds_file")
        self.cov_metric_dropdown = QtWidgets.QComboBox(self.tab)
        self.cov_metric_dropdown.setGeometry(QtCore.QRect(368, 282, 111, 23))
        self.cov_metric_dropdown.setObjectName("cov_metric_dropdown")
        # -- [Add coverage calculation metrics
        self.cov_metric_dropdown.addItems(["Reads/Nt", "RPKM"])
        # -- ]
        self.cov_cutoff_spinbox = QtWidgets.QDoubleSpinBox(self.tab)
        self.cov_cutoff_spinbox.setGeometry(QtCore.QRect(230, 282, 111, 24))
        self.cov_cutoff_spinbox.setObjectName("cov_cutoff_spinbox")
        # -- [Set increment value and default value for coverage cutoff spinbox
        self.cov_cutoff_spinbox.setMaximum(10000)
        self.cov_cutoff_spinbox.setSingleStep(0.05)
        self.cov_cutoff_spinbox.setValue(10)
        # -- ]
        self.listaction_dropdown = QtWidgets.QComboBox(self.tab)
        self.listaction_dropdown.setGeometry(QtCore.QRect(230, 490, 251, 23))
        self.listaction_dropdown.setObjectName("listaction_dropdown")
        # -- [What to do with the gene list
        self.listaction_dropdown.addItems(["Exclusive include", "Exclusive exclude", "Exclude balanced"])
        # -- ]
        self.allowoverlap_entrybox = QtWidgets.QLineEdit(self.tab)
        self.allowoverlap_entrybox.setGeometry(QtCore.QRect(230, 530, 251, 23))
        self.allowoverlap_entrybox.setObjectName("allowoverlap_entrybox")
        self.btn_overlap = QtWidgets.QPushButton(self.tab)
        self.btn_overlap.setGeometry(QtCore.QRect(510, 530, 80, 23))
        self.btn_overlap.setObjectName("btn_overlap")
        # -- [Genes allowed to overlap with each other
        self.btn_overlap.clicked.connect(lambda: self.choose_path_for_fileopen(self.allowoverlap_entrybox))
        # -- ]
        self.operonmember_entrybox = QtWidgets.QLineEdit(self.tab)
        self.operonmember_entrybox.setGeometry(QtCore.QRect(230, 570, 251, 23))
        self.operonmember_entrybox.setObjectName("operonmember_entrybox")
        self.btn_operon = QtWidgets.QPushButton(self.tab)
        self.btn_operon.setGeometry(QtCore.QRect(510, 570, 80, 23))
        self.btn_operon.setObjectName("btn_operon")
        # -- [Select operon member genes
        self.btn_operon.clicked.connect(lambda: self.choose_path_for_fileopen(self.operonmember_entrybox))
        # -- ]
        self.bam_entrybox = QtWidgets.QLineEdit(self.tab)
        self.bam_entrybox.setGeometry(QtCore.QRect(230, 190, 251, 23))
        self.bam_entrybox.setObjectName("bam_entrybox")
        self.lbl_bamfile = QtWidgets.QLabel(self.tab)
        self.lbl_bamfile.setGeometry(QtCore.QRect(30, 193, 171, 16))
        self.lbl_bamfile.setObjectName("lbl_bamfile")
        self.btn_bam = QtWidgets.QPushButton(self.tab)
        self.btn_bam.setGeometry(QtCore.QRect(510, 190, 80, 23))
        self.btn_bam.setObjectName("btn_bam")
        # -- [Select BAM file
        self.btn_bam.clicked.connect(lambda: self.choose_path_for_fileopen(self.bam_entrybox))
        # -- ]
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
        # -- [Leaderless transcript total readcount cutoff
        self.ldtran_sum_cutoff_spinbox.setMaximum(500)
        self.ldtran_sum_cutoff_spinbox.setSingleStep(0.05)
        self.ldtran_sum_cutoff_spinbox.setValue(250)
        # -- ]
        self.lbl_ldtran_sumcut = QtWidgets.QLabel(self.tab_2)
        self.lbl_ldtran_sumcut.setGeometry(QtCore.QRect(48, 186, 181, 16))
        self.lbl_ldtran_sumcut.setObjectName("lbl_ldtran_sumcut")
        self.ldtran_max_cutoff_spinbox = QtWidgets.QDoubleSpinBox(self.tab_2)
        self.ldtran_max_cutoff_spinbox.setGeometry(QtCore.QRect(258, 140, 121, 24))
        self.ldtran_max_cutoff_spinbox.setObjectName("ldtran_max_cutoff_spinbox")
        # --[Leaderless transcript max cutoff]
        self.ldtran_max_cutoff_spinbox.setMaximum(500)
        self.ldtran_max_cutoff_spinbox.setSingleStep(0.05)
        self.ldtran_max_cutoff_spinbox.setValue(40)
        # --]
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
        # --[ Add entries to dropdown
        self.dropdown_readcov_frm_terminal.addItems(["5' end", "3' end"])
        # --]
        self.gtf_transcriptTag_entrybox = QtWidgets.QLineEdit(self.tab_2)
        self.gtf_transcriptTag_entrybox.setGeometry(QtCore.QRect(258, 93, 121, 24))
        self.gtf_transcriptTag_entrybox.setObjectName("gtf_transcriptTag_entrybox")
        self.lbl_transcript_id = QtWidgets.QLabel(self.tab_2)
        self.lbl_transcript_id.setGeometry(QtCore.QRect(48, 91, 191, 21))
        self.lbl_transcript_id.setObjectName("lbl_transcript_id")
        self.ignore_nts_start_spinbox = QtWidgets.QSpinBox(self.tab_2)
        self.ignore_nts_start_spinbox.setGeometry(QtCore.QRect(258, 274, 121, 24))
        self.ignore_nts_start_spinbox.setObjectName("ignore_nts_start_spinbox")
        # --[ignore #nt from 5' end of transcript for calculating gene coverage
        self.ignore_nts_start_spinbox.setMaximum(500)
        self.ignore_nts_start_spinbox.setValue(20)
        # --]
        self.tabWidget.addTab(self.tab_2, "")
        self.buttonBox = QtWidgets.QDialogButtonBox(self.centralwidget)
        self.buttonBox.setGeometry(QtCore.QRect(460, 680, 171, 41))
        self.buttonBox.setStandardButtons(QtWidgets.QDialogButtonBox.Cancel|QtWidgets.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        # --[ OK/Cancel buttons wired up
        self.buttonBox.accepted.connect(lambda: self.clicked_OK())
        self.buttonBox.rejected.connect(lambda: MainWindow.close())
        # --]
        self.btn_export_json = QtWidgets.QPushButton(self.centralwidget)
        self.btn_export_json.setGeometry(QtCore.QRect(36, 688, 131, 23))
        self.btn_export_json.setObjectName("btn_export_json")
        # -- [Export/Save JSON configuration
        self.btn_export_json.clicked.connect(lambda: self.export_config_json())
        # -- ]
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
        # --[ About window
        self.actionAbout.triggered.connect(lambda: self.about_dialog())
        # --]
        self.actionLoad_configuration.setText(_translate("MainWindow", "Load configuration"))
        # -- [ Load JSON config
        self.actionLoad_configuration.triggered.connect(lambda: self.load_json_config())
        self.actionExit_2.setText(_translate("MainWindow", "Exit"))
        # -- [ Exit from menubar
        self.actionExit_2.triggered.connect(lambda: MainWindow.close())

    def run_bard_confirmed(self):
        """
        Confirmation popup before we start
        """
        cnf = QtWidgets.QMessageBox()
        cnf.setIcon(QtWidgets.QMessageBox.Question)
        cnf.setWindowTitle("bard - Query")
        cnf.setText("Do you wish to proceed with the current settings?")
        
        cnf.setStandardButtons(QtWidgets.QMessageBox.Cancel | QtWidgets.QMessageBox.Ok)
        
        response = cnf.exec_()
        
        if response == QtWidgets.QMessageBox.Ok:
            return True
        else:
            return False
        
        
    def about_dialog(self):
        """
        Start bard about dialog window
        """
        abt = QtWidgets.QMessageBox()
        abt.setIcon(QtWidgets.QMessageBox.Information)
        abt.setWindowTitle("About bard")
        abt.setText("This is bard v1.0 (ScriptID=20-12-5-55-42-2020)\t\t\t")
        abnotice = ''' 
        Copyright (C) 2019, 2020  Somdeb Chattopadhyay
        Molecular Genetics Laboratory
        National Institute of Immunology, New Delhi
        
        bard is distributed under the GNU General Public License v3
        '''
        abt.setInformativeText(abnotice)
        abt.exec_()

    def export_config_json(self):
        """
        Triggered by the save configuration button
        """
        file_path_save = self.choose_filepath_for_save()
        
        if file_path_save == "":
            return
        
        self.clicked_OK(exit_when_done=False, confirmation=False)
        save_json(global_config, file_path_save)
        

    def clicked_OK(self, exit_when_done=True, confirmation=True):
        """
        OK button press.
        """
        
        # Ask confirmation first
        if confirmation:
            if self.run_bard_confirmed():
                pass
            else:
                return
        
        # General tab
        global_config["coding_sequence_path"] = self.cds_entrybox.text()
        global_config["coding_sequence_format"] = "fasta" 
        global_config["annotation_file_path"] = self.gtf_entrybox.text()
        global_config["annotation_feature_tag"] = self.gtf_genetag_entrybox.text()
        global_config["annotation_feature_type"] = self.gtf_featuretype_entrybox.text()
        global_config["bam_file_path"] = self.bam_entrybox.text()
        
        if self.offset_rb_5.isChecked():
            global_config["check_offset_from"] = "five_prime"
        if self.offset_rb_3.isChecked():
            global_config["check_offset_from"] = "three_prime"
            
        metric = self.cov_metric_dropdown.currentText()
        if metric == "RPKM":
            global_config["coverage_metric"] = "rpkm"
        if metric == "Reads/Nt":
            global_config["coverage_metric"] = "reads_per_nt"
            
        global_config["coverage_cutoff"] = self.cov_cutoff_spinbox.value()
        
        if self.overlap_ignore_rb_yes.isChecked():
            global_config["will_ignore_overlaps"] = True
        if self.overlap_ignore_rb_no.isChecked():
            global_config["will_ignore_overlaps"] = False
        
        low_range  = int(self.initscan_lower_entrybox.text())
        high_range = int(self.initscan_upper_entrybox.text())
        global_config["peak_scan_range"] = [low_range, high_range]
        
        rdlens = list(
            map(int, self.kmer_entrybox.text().lstrip(",").rstrip(",").strip(" ").split(","))
        )
        global_config["use_readlengths"] = rdlens
        
        global_config["gene_list_file"] = self.genelist_entrybox.text()
        
        action = self.listaction_dropdown.currentText()
        if action == "Exclusive include":
            global_config["gene_list_action"] = "include_only"
        if action == "Exclusive exclude":
            global_config["gene_list_action"] = "exclude_only"
        if action == "Exclude balanced":
            global_config["gene_list_action"] = "exclude_balance"

        global_config["genes_overlap_exception"] = self.allowoverlap_entrybox.text()
        global_config["operon_members_list"] = self.operonmember_entrybox.text()
        
        # Advanced tab
        # We first check if these values have changed from the original
        # preset values first, only then do we modify them
        
        transcript_tag = self.gtf_transcriptTag_entrybox.text()
        ldtran_max = self.ldtran_max_cutoff_spinbox.value()
        ldtran_sum = self.ldtran_sum_cutoff_spinbox.value()
        rdcov_term = self.dropdown_readcov_frm_terminal.currentText()
        rdcov_nts  = self.ignore_nts_start_spinbox.value()
        verbose_filename = self.verbose_filename_checkbox.isChecked()
        
        if transcript_tag != "":
            global_config["transcript_identifier"] = transcript_tag
        
        if ldtran_max != 40:
            global_config["ldtran_max_cutoff"] = ldtran_max
        
        if ldtran_sum != 250:
            global_config["ldtran_sum_cutoff"] = ldtran_sum
            
        if rdcov_nts != 20:
            global_config["readcov_ignore_nts"] = rdcov_nts
            
        if rdcov_term == "3' end":
            global_config["readcov_from_terminal"] = "three_prime"
            
        if verbose_filename:
            global_config["filename_verbosity"] = True
        
        
        if exit_when_done:
            MainWindow.close()

    
    def load_json_config(self):
        """
        Load a saved JSON file for editing
        """
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(MainWindow, 'Load JSON configuration', os.getcwd(),"Configuration files (*.json)")
        
        if fileName == "":
            return
        
        conf = load_json(fileName)
        
        # General tab
        if "coding_sequence_path" in conf:
            self.cds_entrybox.setText(conf["coding_sequence_path"])
            
        if "annotation_file_path" in conf:
            self.gtf_entrybox.setText(conf["annotation_file_path"])
            
        if "annotation_feature_tag" in conf:
            self.gtf_genetag_entrybox.setText(conf["annotation_feature_tag"])
            
        if "annotation_feature_type" in conf:
            self.gtf_featuretype_entrybox.setText(conf["annotation_feature_type"])
            
        if "bam_file_path" in conf:
            self.bam_entrybox.setText(conf["bam_file_path"])
            
        if "check_offset_from" in conf:
            if conf["check_offset_from"] == "five_prime":
                self.offset_rb_5.setChecked(True)
            if conf["check_offset_from"] == "three_prime":
                self.offset_rb_3.setChecked(True)
        
        if "coverage_metric" in conf:
            if conf["coverage_metric"] == "reads_per_nt":
                self.cov_metric_dropdown.setCurrentIndex(0)
            if conf["coverage_metric"] == "rpkm":
                self.cov_metric_dropdown.setCurrentIndex(1)
                
        if "coverage_cutoff" in conf:
            self.cov_cutoff_spinbox.setValue(conf["coverage_cutoff"])
            
        if "will_ignore_overlaps" in conf:
            if conf["will_ignore_overlaps"] == True:
                self.overlap_ignore_rb_yes.setChecked(True)
            if conf["will_ignore_overlaps"] == False:
                self.overlap_ignore_rb_no.setChecked(True)
        
        if "use_readlengths" in conf:
            rdlen_str = ",".join(list(map(str, conf["use_readlengths"])))
            self.kmer_entrybox.setText(rdlen_str)
            
        if "peak_scan_range" in conf:
            lowerlimit = str(conf["peak_scan_range"][0])
            upperlimit = str(conf["peak_scan_range"][1])
            
            self.initscan_lower_entrybox.setText(lowerlimit)
            self.initscan_upper_entrybox.setText(upperlimit)
            
        if "gene_list_file" in conf:
            genelist = conf["gene_list_file"]
            self.genelist_entrybox.setText(genelist)
            
        if "gene_list_action" in conf:
            action = conf["gene_list_action"]
            
            if action == "include_only":
                self.listaction_dropdown.setCurrentIndex(0)
            if action == "exclude_only":
                self.listaction_dropdown.setCurrentIndex(1)
            if action == "exclude_balance":
                self.listaction_dropdown.setCurrentIndex(2)
                
        if "genes_overlap_exception" in conf:
            overlap_allow = conf["genes_overlap_exception"]
            self.allowoverlap_entrybox.setText(overlap_allow)
            
        if "operon_members_list" in conf:
            omlist = conf["operon_members_list"]
            self.operonmember_entrybox.setText(omlist)
            
        # Advanced tab
        if "transcript_identifier" in conf:
            tid = conf["transcript_identifier"]
            self.gtf_transcriptTag_entrybox.setText(tid)
            
        if "ldtran_max_cutoff" in conf:
            lmcutoff = conf["ldtran_max_cutoff"]
            self.ldtran_max_cutoff_spinbox.setValue(lmcutoff)
            
        if "ldtran_sum_cutoff" in conf:
            lscutoff = conf["ldtran_sum_cutoff"]
            self.ldtran_sum_cutoff_spinbox.setValue(lscutoff)
            
        if "readcov_ignore_nts" in conf:
            rdcov_ign = conf["readcov_ignore_nts"]
            self.ignore_nts_start_spinbox.setValue(rdcov_ign)
            
        if "readcov_from_terminal" in conf:
            rcterm = conf["readcov_from_terminal"]
            
            if rcterm == "five_prime":
                self.dropdown_readcov_frm_terminal.setCurrentIndex(0)
            if rcterm == "three_prime":
                self.dropdown_readcov_frm_terminal.setCurrentIndex(1)
                
        if "filename_verbosity" in conf:
            verbose = conf["filename_verbosity"]
            self.verbose_filename_checkbox.setChecked(verbose)
            
    
    def choose_path_for_fileopen(self, fobject):
        """
        we set the file path in the entrybox,
        we then retrieve the path from the entrybox itself
        connected to the entrybox corresponding button
        """
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(MainWindow, 'Open file', os.getcwd())
        
        if fileName == "":
            return
        
        fobject.setText(fileName)
        
    def choose_filepath_for_save(self):
        """
        GUI dialog to choose location to save a file
        (JSON configuration)
        """
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(MainWindow, 'Save JSON configuration', os.getcwd(),"Configuration files (*.json)")
        #self.file_paths[filetype] = fileName
        return(fileName)

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())

