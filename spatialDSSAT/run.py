import pickle
import numpy as np

import os
import shutil
import subprocess
import re
import string
from itertools import product
from tqdm import tqdm
import tempfile 

from utils import *

from datetime import datetime
import random

import warnings
warnings.filterwarnings(action='ignore')

TMP =  tempfile.gettempdir()
WTH_IDS = ["".join(i) for i in product(string.ascii_uppercase, string.ascii_uppercase)]

# RUN_PATH = f"/tmp/dssatrun"
# if not os.path.exists(RUN_PATH):
#     os.mkdir(RUN_PATH)

HOME = "~"
# DSSAT_HOME = "/home/dquintero/DSSAT48/"
# DSSAT_BIN = "/home/dquintero/dssat-csm-os/build/bin/dscsm048"
COUNTRY = "Kenya"


class GSRun():
    '''
    A class to handle the DSSAT execution environment. The DSSAT execution environment
    is the folder where all the inputs and outputs will be located. When this classes
    is instantiated then a tmp folder is created to run the model.

    It runs DSSAT in the Spatial mode. This mode basically runs individual experiments
    for the different locations. See: https://dssat.net/wp-content/uploads/2011/10/DSSAT-vol4.pdf
    for further details.
    '''
    def __init__(self, dssat_home=None, dssat_bin=None):
        '''
        Instantiate the GSRun class

        Arguments
        ----------
        dssat_home: str
            Path to the folder with the files DSSAT needs. Those files include 
            .CDE, .CTR, .L48 files, and Genotype, Pest, and StandardData. If not
            provided it is taken from DSSATTools path.
        dssat_bin: str
            Path to the DSSAT executable file. If not provided it is taken from
            the DSSATTools path.
        '''
        self.__DSSAT_BIN = dssat_bin
        self.__DSSAT_HOME = dssat_home

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
        f'! Command Line : {self.__DSSAT_BIN} S DSSBatch.v48\n' + \
        f"! Crop         : Spatial\n" + \
        f"! Experiment   : EXPEFILE.MZX\n" + \
        f"! ExpNo        : 1\n" + \
        f'! Debug        : {self.__DSSAT_BIN} " S DSSBatch.v48"\n' + \
        "!\n" + \
        "@FILEX                                                                                        TRTNO     RP     SQ     OP     CO\n"
        
    def add_treatment(self, soil:tuple, weather:tuple, nitrogen:int, 
                      planting:datetime, cultivar:str):
        """
        It adds a treatment for the current run.
        """
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
            f"{os.path.join(self.RUN_PATH, 'EXPEFILE.MZX'):<96} {n:-2}      1      0      0      0\n"
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
            soil = self.soil[soil_n-1]
            soil_profile = f"IB000000{soil_n:-02}"
            soil_filename = f"{soil[0]:07.2f}_{soil[1]:07.2f}".replace(".", "p")+".SOL"
            with open(f"{HOME}/data/soil_data/dssat_soils/{soil_filename}", "r") as f:
                soil_lines = f.readlines()

            fieldCap = []
            for layer_line in soil_lines[6:]:
                if len(layer_line) > 2:
                    fieldCap.append((layer_line[:6], layer_line[19:24]))
            fieldCap = np.array(fieldCap).astype(float)
            fieldCap = np.hstack([fieldCap, np.ones(fieldCap.shape)*.01]) # Default value for initial NO3 and NH4

            "@C   PCR ICDAT  ICRT  ICND  ICRN  ICRE  ICWD ICRES ICREN ICREP ICRIP ICRID ICNAME"
            " 4    MZ 81089   200   -99     1     1   -99   -99   -99   -99   -99   -99 -99"
            "@C  ICBL  SH2O  SNH4  SNO3"
            " 4    30  .107    .2    .7"

            # management.initial_conditions["table"] = TabularSubsection(pd.DataFrame(
            #     fieldCap,
            #     columns=['ICBL', 'SH2O', 'SNH4', 'SNO3'])
            # )
            # management.field['SLDP'] = fieldCap[:,0].max()
            # management.field["...........XCRD"] = soil[0]
            # management.field["...........YCRD"] = soil[1]
    
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

    def _options_build(self):
        self.gsx_str += \
        "*SIMULATION CONTROLS\n" + \
        "@N GENERAL     NYERS NREPS START SDATE RSEED SNAME.................... SMODEL\n" + \
        " 1 GE              1     1     P {0:<5}  2150 N SPATIAL ANALYSES TEST\n".format(
            self.start_date.strftime("%y%j")
        ) + \
        "@N OPTIONS     WATER NITRO SYMBI PHOSP POTAS DISES  CHEM  TILL   CO2\n" + \
        " 1 OP              Y     Y     N     N     N     N     N     N     D\n" + \
        "@N METHODS     WTHER INCON LIGHT EVAPO INFIL PHOTO HYDRO NSWIT MESOM MESEV MESOL\n" + \
        " 1 ME              M     M     E     R     S     C     R     1     P     R     2\n" + \
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

    def run(self):
        self._treatment_build()
        self._cultivar_build()
        self._field_build()
        self._planting_build()
        self._fertilizer_build()
        self._options_build()
        # self._ic_build()

        for n, soil in enumerate(self.soil, 1):
            soil_filename = f"{soil[0]:07.2f}_{soil[1]:07.2f}".replace(".", "p")+".SOL"
            with open(f"{HOME}/data/soil_data/dssat_soils/{soil_filename}", "r") as f:
                soil_lines = f.readlines()
            soil_lines[0] = f"*IB000000{n:-02}" + soil_lines[0][11:]
            self.sol_str += "".join(soil_lines)
        with open(f"{self.RUN_PATH}/SOIL.SOL", 'w') as f:
            f.write(self.sol_str)

        for n, weather in enumerate(self.weather, 1):
            wth_id = WTH_IDS[n-1]
            for year in range(self.start_date.year-1, self.start_date.year+1):
                wthfile_path = f"{self.RUN_PATH}/SE{wth_id}{str(year)[2:]}01.WTH"
                if os.path.exists(wthfile_path):
                    os.remove(wthfile_path)
                wthfile = f"{weather[0]:07.2f}_{weather[1]:07.2f}_{str(year)[2:]}01".replace(".", "p")+".WTH"    
                assert os.path.exists(f"{HOME}/data/weather_data/dssat_files/{COUNTRY}/{wthfile}")
                os.symlink(
                    f"{HOME}/data/weather_data/dssat_files/{COUNTRY}/{wthfile}",
                    wthfile_path
                )

        with open(f"{self.RUN_PATH}/EXPEFILE.MZX", 'w') as f:
            f.write(self.gsx_str)

        with open(f"{self.RUN_PATH}/DSSBatch.v48", 'w') as f:
            f.write(self.batch_str)

        exc_args = [f"{self.__DSSAT_BIN}", 'S', "DSSBatch.v48"]
        excinfo = subprocess.run(exc_args, 
            cwd=self.RUN_PATH, capture_output=True, text=True,
            env={"DSSAT_HOME": self.__DSSAT_HOME, }
        )
        out = re.sub(r'(\n{2,})|(\n$)', "", excinfo.stdout)
        out = re.sub(r'((RUN).+\n.+(t/ha)\n)', "", out)
        out = re.sub(f'\n.+(Crop).+\n', "\n", out)
        df = pd.DataFrame(
            [line.split() for line in out.split("\n")],
            columns="RUN  CR  TRT FLO MAT TOPWT HARWT  RAIN  TIRR   CET  PESW  TNUP  TNLF   TSON TSOC".split()
        )
        shutil.rmtree(self.RUN_PATH)
        return df
    

if __name__ == "__main__":
    planting_window = pd.date_range("2014-03-20", "2014-05-20") # Long rains planting

    with open(f"data/pixels_counties.pkl", "rb") as f:
        pixels = pickle.load(f)
    cultivars = ["990001", "990002", "990003"]
    PRODUCTIVE_COUNTIES = ['Trans Nzoia', 'Uasin Gishu', 'Bungoma', 'Kakamega', 'Narok', 'Nakuru',
        'Nandi', 'Migori', 'Kisii', 'Kericho', 'Machakos', 'Elgeyo-Marakwet',
        'Siaya', 'Homa Bay', 'Meru', 'Makueni', 'Kwale', 'West Pokot', 'Kilifi',
        'Baringo', 'Nyamira', 'Busia', 'Muranga', 'Bomet', 'Kisumu']
    # PRODUCTIVE_COUNTIES = ["Uasin Gishu"]
    N = 99
    df_runs = pd.DataFrame()
    for year in range(2014, 2020):
        for province, pix in tqdm(pixels.items()):
            if province not in PRODUCTIVE_COUNTIES:
                continue
            soil_list = [tuple(pix[np.random.randint(len(pix))]) for _ in range(N)]
            weather_list = [tuple(pix[np.random.randint(len(pix))]) for _ in range(N)]
            nitrogen_list = [((0, 0), )]*N
            planting_list = [
                datetime(year, planting_window[i].month, planting_window[i].day)
                for i in np.random.randint(len(planting_window), size=N)
            ]
            cultivars_list = [cultivars[np.random.randint(0, 3)] for _ in range(N)]

            gsRun = GSRun()
            for treat in zip(soil_list, weather_list, nitrogen_list, planting_list, cultivars_list):    
                gsRun.add_treatment(*treat)
            out = gsRun.run()
            for file in os.listdir(gsRun.RUN_PATH):
                if file[-3:] in ("WTH", "SOL"):
                    os.remove(os.path.join(gsRun.RUN_PATH, file))
            out["soil"] = soil_list 
            out["weather"] = weather_list 
            out["nitro"] = nitrogen_list 
            out["planting"] = planting_list 
            out["cultivars"] = cultivars_list
            out["year"] = year
            out["province"] = province

            df_runs = pd.concat([df_runs, out], ignore_index=True)

        # df_runs.to_csv("results_tmp.csv")
    df_runs.to_csv("results.csv")
    # df_runs.to_csv("results_debug.csv")