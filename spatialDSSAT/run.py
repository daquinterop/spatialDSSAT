import numpy as np

import os
import shutil
import subprocess
import re
import string
from itertools import product
import tempfile 

from spatialDSSAT.utils import *
from DSSATTools import __file__ as dsssattools_module_path
from DSSATTools import VERSION
from DSSATTools.crop import CROP_CODES, CROPS_MODULES

from datetime import datetime, timedelta
import random

import warnings
warnings.filterwarnings(action='ignore')

TMP =  tempfile.gettempdir()
# String squence to name weather files (AA, AB, AC, ..., XZ, YZ, ZZ)
WTH_IDS = ["".join(i) for i in product(string.ascii_uppercase, string.ascii_uppercase)]

DSSAT_STATIC = os.path.join(os.path.dirname(dsssattools_module_path), "static")
DSSAT_BIN = os.path.join(
    os.path.dirname(dsssattools_module_path), "static", "bin", "dscsm048"
)
DSSAT_HOME = os.path.join(TMP, f"DSSAT{VERSION}/")
BIN_NAME = f"dscsm{VERSION}"
CONFILE = f'DSSATPRO.L{VERSION[1:]}'

# Creates a folder with DSSAT files. This is done to avoid long path names that 
# exceed the defined lenght for path variables in DSSAT.
if not os.path.exists(DSSAT_HOME):
    os.mkdir(DSSAT_HOME)
for file in os.listdir(DSSAT_STATIC):
    file_link = os.path.join(DSSAT_HOME, file)
    if os.path.exists(file_link):
        os.remove(file_link)
    os.symlink(os.path.join(DSSAT_STATIC, file), file_link)

DEFAULT_SIMULATION_OPTIONS = {
    "WATER": "Y", "NITRO": "Y", "SYMBI": "N", "PHOSP": "N", "POTAS": "N",
    "CO2": "D",
    "LIGHT": "E", "EVAPO": "R", "INFIL": "S", "PHOTO": "C",  "MESOM": "P", 
    "MESEV": "R", "MESOL": "2",
}

def write_control_file(run_path, crop_name):
    """
    Writes DSSAT control file in the specified directory.

    Arguments
    ----------
    run_path: str
        Path to the directory where the model will run
    crop_code: str
        Crop name
    """
    crop_code = CROP_CODES[crop_name]
    smodel = CROPS_MODULES[crop_name]
    wth_path = "Weather"
    with open(os.path.join(run_path, CONFILE), 'w') as f:
        f.write(f'WED    {wth_path}\n')
        if crop_code in ["WH", "BA"]:
            f.write(f'M{crop_code}    {run_path} dscsm048 CSCER{VERSION}\n')
        else:
            f.write(f'M{crop_code}    {DSSAT_HOME} dscsm048 {smodel}{VERSION}\n')
        f.write(f'CRD    {os.path.join(DSSAT_HOME, "Genotype")}\n')
        f.write(f'PSD    {os.path.join(DSSAT_HOME, "Pest")}\n')
        f.write(f'SLD    {os.path.join(DSSAT_HOME, "Soil")}\n')
        f.write(f'STD    {os.path.join(DSSAT_HOME, "StandardData")}\n')


class GSRun():
    '''
    A class to handle the DSSAT execution environment. The DSSAT execution environment
    is the folder where all the inputs and outputs will be located. When this classes
    is instantiated then a tmp folder is created to run the model.

    It runs DSSAT in the Spatial mode. This mode basically runs individual experiments
    for the different locations. See: https://dssat.net/wp-content/uploads/2011/10/DSSAT-vol4.pdf
    for further details.
    '''
    def __init__(self, crop_name="Maize"):
        '''
        Instantiate the GSRun class

        Arguments
        ----------
        crop_name: str
            Name of the crop. Default value is Maize
        '''
        if not hasattr(self, "treatments"):
            self._crop_name = crop_name.title()
            assert self._crop_name in CROP_CODES.keys(), \
                f'{self._crop_name} is not a valid crop'

            self.RUN_PATH = os.path.join(
                TMP, "dssatrun"+''.join(random.choices(string.ascii_uppercase + string.digits, k=8))
            )
            if not os.path.exists(self.RUN_PATH):
                os.mkdir(self.RUN_PATH)      
        
        self.weather = []
        self.soil = []
        self.field = []
        self.planting = []
        self.nitrogen = []
        self.cultivar = []
        self.treatments = []
        self.start_date = datetime(9999, 9, 9)
        self.sol_str = ""
        self.gsx_str = ""
        self.batch_str = ""
        
    def add_treatment(self, soil:str, weather:str, nitrogen:list, 
                      planting:datetime, cultivar:str):
        """
        It adds a treatment (location) for the current run. It adds one treatment at a time.

        Arguments
        ----------
        soil: str 
            Path to a .SOL file. That file must contain only one soil profile.
        weather: str
            Path to a .WTH file. 
        nitrogen: list of tuples [(int, float), ...]
            A list of tuples where each tuple is one Nitrogen application. The 
            tuple contains two values, where the first one indicates the application
            time (days after planting) and the second one the nitrogen rate (kg/ha).
        planting: datetime
            Planting date
        cultivar: str
            Cultivar code. It must be available on the DSSAT list of cultivars.
        """

        assert len(self.treatments) < 99, "99 is the maximum number of treatments"

        if weather not in self.weather:
            self.weather.append(weather)
        weather_n = self.weather.index(weather) + 1

        if soil not in self.soil:
            self.soil.append(soil)
        soil_n = self.soil.index(soil) + 1

        field = (weather_n, soil_n)
        if field not in self.field:
            self.field.append(field)
        field_n = self.field.index(field) + 1

        if nitrogen not in self.nitrogen:
            self.nitrogen.append(nitrogen)
        nitrogen_n = self.nitrogen.index(nitrogen) + 1

        if cultivar not in self.cultivar:
            self.cultivar.append(cultivar)
        cultivar_n = self.cultivar.index(cultivar) + 1

        if planting not in self.planting:
            self.planting.append(planting)
        planting_n = self.planting.index(planting) + 1
        self.start_date = min(self.start_date, planting)

        self.treatments.append((cultivar_n, field_n, planting_n, nitrogen_n))

    def _header_build(self):
        self.sol_str = "*SOILS: the upper layer of earth in which plants grow\n\n"
        self.gsx_str = \
        "*EXP.DETAILS: FLSC8101GS SPATIAL ANALYSES TEST CASE; FLORENCE, SOUTH CAROLINA\n\n" +\
        "*GENERAL\n" +\
        "@PEOPLE\n" +\
        "Diego\n" +\
        "@ADDRESS\n" +\
        "Huntsville, Alabama\n" +\
        "@SITE\n" +\
        "Site\n" +\
        "@ PAREA  PRNO  PLEN  PLDR  PLSP  PLAY HAREA  HRNO  HLEN  HARM.........\n" +\
        "    -99   -99   -99   -99   -99   -99   -99   -99   -99   -99\n\n"
        self.batch_str = \
        "$BATCH(SPATIAL)\n" + \
        "!\n" + \
        f"! Directory    : {self.RUN_PATH}\n" + \
        f'! Command Line : {BIN_NAME} S DSSBatch.v48\n' + \
        f"! Crop         : Spatial\n" + \
        f"! Experiment   : EXPEFILE.{CROP_CODES[self._crop_name]}X\n" + \
        f"! ExpNo        : 1\n" + \
        f'! Debug        : {BIN_NAME} " S DSSBatch.v48"\n' + \
        "!\n" + \
        "@FILEX                                                                                        TRTNO     RP     SQ     OP     CO\n"

    def _treatment_build(self):
        self.gsx_str += \
        "*TREATMENTS                        -------------FACTOR LEVELS------------\n" +\
        "@N R O C TNAME.................... CU FL SA IC MP MI MF MR MC MT ME MH SM\n" 
        for n, treat in enumerate(self.treatments, 1):
            cultivar_n, field_n, planting_n, nitrogen_n = treat
            self.gsx_str += \
            "{:-2} 1 0 0 SAMPLE {:<18} {:-2} {:-2}  0 {:-2} {:-2}  0 {:-2}  0  0  0  0  0  1\n".format(
                n, n, cultivar_n, field_n, 0, planting_n, nitrogen_n
            ) # IC is 0, then default options are used (Field capacity, zero nitrogen)
            self.batch_str += \
            f"{os.path.join(self.RUN_PATH, f'EXPEFILE.{CROP_CODES[self._crop_name]}X'):<96} {n:-2}      1      0      0      0\n"
        self.gsx_str += "\n" 

    def _cultivar_build(self):
        self.gsx_str += \
        "*CULTIVARS\n" +\
        "@C CR INGENO CNAME\n"
        for n, cul in enumerate(self.cultivar, 1):
            self.gsx_str += \
            "{:-2} MZ {:>6} {:<8}\n".format(n, cul, cul)
        self.gsx_str += "\n" 

    def _field_build(self):
        yr = str(self.start_date.year)[2:]
        self.gsx_str += \
        "*FIELDS\n" + \
        "@L ID_FIELD WSTA....  FLSA  FLOB  FLDT  FLDD  FLDS  FLST SLTX  SLDP  ID_SOIL    FLNAME\n" 
        for n, (weather_n, soil_n) in enumerate(self.field, 1):
            wth_id = WTH_IDS[weather_n - 1]
            self.gsx_str += \
            "{0:-2} SEFL00{0:02} SE{1}{2}01   -99     0 DR000     0     0 00000 -99    200  IB000000{3:-02} -99\n".format(n, wth_id, yr, soil_n)
        self.gsx_str += "\n" 
        
    def _ic_build(self):
        """
        Build the initial conditions section. Creates the soil file.
        """
        for n, (_, field_n, _, _) in enumerate(self.treatments, 1):
            soil_n = self.field[field_n-1][1]
            soil_path = self.soil[soil_n-1]
            with open(soil_path, "r") as f:
                soil_lines = f.readlines()

            fieldCap = []
            for layer_line in soil_lines[6:]:
                if len(layer_line) > 2:
                    fieldCap.append((layer_line[:6], layer_line[19:24]))
            fieldCap = np.array(fieldCap).astype(float)
            fieldCap = np.hstack([fieldCap, np.ones(fieldCap.shape)*.01]) # Default value for initial NO3 and NH4

    
    def _planting_build(self):
        self.gsx_str += \
        "*PLANTING DETAILS\n" + \
        "@P PDATE EDATE  PPOP  PPOE  PLME  PLDS  PLRS  PLRD  PLDP  PLWT  PAGE  PENV  PLPH  SPRL                        PLNAME\n"
        for n, planting in enumerate(self.planting, 1):
            self.gsx_str += \
            "{0:-2} {1:<5}   -99   4.0   4.0     S     R    90     0     4   -99   -99   -99   -99   -99                        -99\n". format(
                n, planting.strftime("%y%j")
            ) 
        self.gsx_str += "\n" 

    def _fertilizer_build(self):
        self.gsx_str += \
        "*FERTILIZERS (INORGANIC)\n" + \
        "@F FDATE  FMCD  FACD  FDEP  FAMN  FAMP  FAMK  FAMC  FAMO  FOCD FERNAME\n" 
        for nitrogen_n, nitrogen in enumerate(self.nitrogen, 1):
            for dap, rate in nitrogen:
                self.gsx_str += \
                "{0:-2} {1:<5} FE005   -99     5 {2:5.1f}   -99   -99   -99   -99   -99 -99\n".format(
                    nitrogen_n, dap, rate
                )
        self.gsx_str += "\n" 

    def _options_build(self, sim_controls):
        options = [
            sim_controls.get(key, DEFAULT_SIMULATION_OPTIONS[key]) 
            for key in DEFAULT_SIMULATION_OPTIONS
        ]
        self.gsx_str += \
        "*SIMULATION CONTROLS\n" + \
        "@N GENERAL     NYERS NREPS START SDATE RSEED SNAME.................... SMODEL\n" + \
        " 1 GE              1     1     P {0:<5}  2150 N SPATIAL ANALYSES TEST\n".format(
            self.start_date.strftime("%y%j")
        ) + \
        "@N OPTIONS     WATER NITRO SYMBI PHOSP POTAS DISES  CHEM  TILL   CO2\n" + \
        " 1 OP              {0}     {1}     {2}     {3}     {4}     N     N     N     {5}\n".format(*options) + \
        "@N METHODS     WTHER INCON LIGHT EVAPO INFIL PHOTO HYDRO MESOM MESEV MESOL\n" + \
        " 1 ME              M     M     {6}     {7}     {8}     {9}     R     {10}     {11}     {12}\n".format(*options) + \
        "@N MANAGEMENT  PLANT IRRIG FERTI RESID HARVS\n" + \
        " 1 MA              R     N     D     N     M\n" + \
        "@N OUTPUTS     FNAME OVVEW SUMRY FROPT GROUT CAOUT WAOUT NIOUT MIOUT DIOUT VBOSE CHOUT OPOUT FMOPT\n" + \
        " 1 OU              N     N     N     1     N     N     N     N     N     N     Y     N     N     A\n\n" + \
        "@  AUTOMATIC MANAGEMENT\n" + \
        "@N PLANTING    PFRST PLAST PH2OL PH2OU PH2OD PSTMX PSTMN\n" + \
        " 1 PL          98155 98200    40   100    30    40    10\n" + \
        "@N IRRIGATION  IMDEP ITHRL ITHRU IROFF IMETH IRAMT IREFF\n" + \
        " 1 IR             30    50   100 GS000 IR001    10     1\n" + \
        "@N NITROGEN    NMDEP NMTHR NAMNT NCODE NAOFF\n" + \
        " 1 NI             30    50    25 FE001 GS000\n" + \
        "@N RESIDUES    RIPCN RTIME RIDEP\n" + \
        " 1 RE            100     1    20\n" + \
        "@N HARVEST     HFRST HLAST HPCNP HPCNR\n" + \
        " 1 HA              0 81365   100     0\n"

    def run(self, **kwargs) -> pd.DataFrame:
        """
        Run DSSAT in spatial mode. It returns a dataframe with the simulation
        reults. No arguments are required, some are optional tough.

        Arguments
        ----------
        start_date: datetime
            Start date for all treatments. If not provided the earliest planting
            date is taken.
        latest_date: datetime
            Latest date the simulation is expected to end. This parameter is used 
            to avoid searching weather files that might not be available.
        sim_controls: dict
            Simulation control definitions. A dict defining some of the simulation
            options. An example can be found in spatialDSSAT.run.DEFAULT_SIMULATION_OPTIONS
        """
        assert len(self.treatments) > 0, \
            "No treatments have been added. Use the add_treatment to add treatments" 
        self.start_date = kwargs.get("start_date", self.start_date + timedelta(days=150))
        latest_date = kwargs.get("latest_date", self.start_date + timedelta(days=150))
        sim_controls = kwargs.get("sim_controls", DEFAULT_SIMULATION_OPTIONS)
        
        self._header_build()
        self._treatment_build()
        self._cultivar_build()
        self._field_build()
        self._planting_build()
        self._fertilizer_build()
        self._options_build(sim_controls)
        
        write_control_file(self.RUN_PATH, self._crop_name)

        for n, soil_path in enumerate(self.soil, 1):
            with open(soil_path, "r") as f:
                soil_lines = f.readlines()
            soil_lines[0] = f"*IB000000{n:-02}" + soil_lines[0][11:]
            self.sol_str += "".join(soil_lines)
        with open(f"{self.RUN_PATH}/SOIL.SOL", 'w') as f:
            f.write(self.sol_str)

        for n, weather_path in enumerate(self.weather, 1):
            wth_id = WTH_IDS[n-1]
            for year in range(self.start_date.year, latest_date.year+1):
                wthfile_path = f"{self.RUN_PATH}/SE{wth_id}{str(year)[2:]}01.WTH"
                if os.path.exists(wthfile_path):
                    os.remove(wthfile_path)
                assert os.path.exists(weather_path)
                os.symlink(weather_path, wthfile_path)

        with open(f"{self.RUN_PATH}/EXPEFILE.MZX", 'w') as f:
            f.write(self.gsx_str)

        with open(f"{self.RUN_PATH}/DSSBatch.v48", 'w') as f:
            f.write(self.batch_str)

        exc_args = [f"{DSSAT_BIN}", 'S', "DSSBatch.v48"]
        excinfo = subprocess.run(exc_args, 
            cwd=self.RUN_PATH, capture_output=True, text=True,
            env={"DSSAT_HOME": DSSAT_HOME, BIN_NAME: DSSAT_BIN}
        )
        out = re.sub(r'(\n{2,})|(\n$)', "", excinfo.stdout)
        out = re.sub(r'((RUN).+\n.+(t/ha)\n)', "", out)
        out = re.sub(f'\n.+(Crop).+\n', "\n", out)
        df = pd.DataFrame(
            [line.split() for line in out.split("\n")],
            columns="RUN  CR  TRT FLO MAT TOPWT HARWT  RAIN  TIRR   CET  PESW  TNUP  TNLF   TSON TSOC".split()
        )
        # shutil.rmtree(self.RUN_PATH)
        return df
    
    def clear(self):
        """
        Clear treatments to define new ones
        """
        self.__init__(crop_name=self._crop_name)