import datetime
import csv
import collections
import numpy
from osgeo import osr
from osgeo import gdal
import logging
import sys
from ecoshard import taskgraph
import os
import geopandas as gpd
from ecoshard import geoprocessing

logging.basicConfig(
    level=logging.DEBUG,
    stream=sys.stdout,
    format=(
        '%(asctime)s (%(relativeCreated)d) %(levelname)s %(name)s'
        ' [%(pathname)s.%(funcName)s:%(lineno)d] %(message)s'))
LOGGER = logging.getLogger(__name__)
logging.getLogger('matplotlib.font_manager').setLevel(logging.ERROR)
logging.getLogger('PIL').setLevel(logging.ERROR)
logging.getLogger('ecoshard.taskgraph').setLevel(logging.INFO)
logging.getLogger('fiona').setLevel(logging.WARN)

#DEM_PATH = r"D:/repositories/downstream-beneficiaries/workspace/global_dem_3s_md5_22d0c3809af491fa09d03002bdf09748/global_dem_3s"

PERCENTILES_LIST = []

# report the area in the clipped raster that has values >= to these values
THRESHOLD_AREA_LIST = []

BASE_RASTER_LOOKUP = {
    #'lspop2019': r"D:\repositories\roadmap2030\data/pop_rasters/lspop2019_compressed_md5_d0bf03bd0a2378196327bbe6e898b70c.tif",
    #'floodpop': r"D:\repositories\roadmap2030\data/pop_rasters/floodplains_masked_pop_30s_md5_c027686bb9a9a36bdababbe8af35d696.tif",
    #'lspop2023': r"D:\repositories\roadmap2030\data\pop_rasters\landscan-global-2023.tif",
    'eii': r"D:\repositories\data_platform\Nature\eii_padj_v5140524_epsg_3395.tif",
    #'sed_export_change': r"D:\repositories\roadmap2030\data\ndv_0.0_sed_export_marineESA_2020-1992_change_md5_0ab0cf.tif",
    ##'cv_habitat_change': r"D:\repositories\roadmap2030\data\ndv_0.0_cv_habitat_value_marESA2020-1992_change_md5_1643a7.tif",
    #'n_export_change': r"D:\repositories\roadmap2030\data\ndv_0.0_n_export_marineESA_2020-1992_change_val_md5_18a2b3.tif",
    #'realized_pollination_on_ag_change': r"D:\repositories\roadmap2030\data\ndv_0.0_realized_pollination_on_ag_marESA_2020-1992_fullchange_md5_8e63e2.tif",
    #'sed_deposition_change': r"D:\repositories\roadmap2030\data\ndv_0.0_sed_deposition_marineESA_2020-1992_change_md5_d23c49.tif",
    ##'coastal_reference': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\coastal_risk_Sc3v1_habitat_value_md5_e889c2dbc5783fc4c782fbd3b473d7de.tif",
    ##'coastal_change': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\coastal_risk_tnc_esa2020_change_esa1992_md5_ea900e.tif",
    ##'coastal_2020': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\coastal_risk_tnc_esa2020_value_md5_f9f644.tif",
    #'nitrogen_reference_full': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\n_export_sc3v2pnvall_compressed_md5_09bc65fe1cd54b518cde859f57513d8c.tif",
    #'nitrogen_reference': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\n_export_sc3v1pnvnoag_compressed_md5_bd5a856e0c1f76b2e8898f533ec20659.tif",
#    'nitrogen_change': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\n_export_tnc_2020-1992_change_val_md5_18a2b3.tif",
#    'nitrogen_2020': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\n_export_tnc_esa2020_compressed_md5_1d3c17.tif",
#    'pollination_change': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\pollination_on_ag_marESA_2020-1992_fullchange_md5_8e63e2.tif",
#    #'pollination_reference': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\pollination_ppl_fed_on_ag_10s_Sc3v1_PNVnoag.tif",
#    'pollination_2020v2': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\pollination_ppl_fed_on_ag_10s_tnc_esa2020ag_compressed_md5_8b5ee8.tif",
#    'pollination_2020v1': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\polllination_on_ag_ESA2020mar_md5_da610a.tif",
#    #'sediment_reference_full': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\sed_export_pnv_compressed_md5_a1faed.tif",
#    #'sediment_reference': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\sed_export_sc3v1pnvnoag_compressed_md5_2783ee50e908a763622d3167669b60bc.tif",
#    'sediment_change': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\sed_export_tnc_ESA_2020-1992_change_md5_0ab0cf.tif",
#    'sediment_2020': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\sed_export_tnc_ESA_2020_compressed_md5_a988c0.tif",
#    #'n_retention_2020': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\n_retention_esamod2_compressed_md5_30d56daec1140d031aa62a2bd6fe1f63.tif",
#    #'n_retention_reference': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\n_retention_sc3v1pnvnoag_compressed_md5_ffb5af20f07c64deb67fcd2e0ffd628d.tif",
#    #'sed_retention_2020': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\sed_retention_esamod2_compressed_md5_c7a77e50feaea7a5dc7322cd63f0f429.tif",
#    #'sed_retention_ref': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\sed_retention_sc3v1pnvnoag_compressed_md5_f6def2b90f231703e813dea293e67fd2.tif",
#    #'sed_retention_ref_full': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\sed_retention_sc3v2pnvall_compressed_md5_f3c11ea7d473237be2b1dadaa9efb172.tif",
#    #'sed_deposition_ref': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\sed_deposition_sc3v1pnvnoag_compressed_md5_d6eeb1717eff44d89f19a57dceab6328.tif",
#    #'sed_deposition_ref_full': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\sed_deposition_sc3v2pnvall_compressed_md5_4c7847e31c5f7afc9ef39e3abec1912d.tif",
#    #'sed_deposition_2020': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\sed_deposition_esamod2_compressed_md5_ff134776cd7d9d69dc5e2fe14b53474c.tif",
#    #'flood_2020': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\flood_nathab_md5_eb8fd58621e00c6aeb80f4483da1b35c.tif",
#    #'poll_suff_2020': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\poll_suff_ag_coverage_prop_10s_ESACCI-LC-L4-LCCS-Map-300m-P1Y-2020-v2.1.1_md5_2ed6285e6f8ec1e7e0b75309cc6d6f9f.tif",
#    #'carbon_reference_full': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\restoration_limited_md5_372bdfd9ffaf810b5f68ddeb4704f48f_forest_projected_full_forest_edge_result.tif",
#    #'carbon_reference': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\restoration_limited_md5_372bdfd9ffaf810b5f68ddeb4704f48f_forest_projected_no_forest_edge_result.tif",
#    'carbon_2000': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\fc_stack_hansen_forest_cover2000_compressed_std_forest_edge_result.tif",
#    'carbon_2020': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\fc_stack_hansen_forest_cover2020_compressed_std_forest_edge_result.tif",
    #'abaca': r"D:\repositories\roadmap2030\data\crop_maps\abaca_harea.tif",
    #'agave': r"D:\repositories\roadmap2030\data\crop_maps\agave_harea.tif",
    #'alfalfa': r"D:\repositories\roadmap2030\data\crop_maps\alfalfa_harea.tif",
    #'almond': r"D:\repositories\roadmap2030\data\crop_maps\almond_harea.tif",
    #'aniseetc': r"D:\repositories\roadmap2030\data\crop_maps\aniseetc_harea.tif",
    #'apple': r"D:\repositories\roadmap2030\data\crop_maps\apple_harea.tif",
    #'apricot': r"D:\repositories\roadmap2030\data\crop_maps\apricot_harea.tif",
    #'areca': r"D:\repositories\roadmap2030\data\crop_maps\areca_harea.tif",
    #'artichoke': r"D:\repositories\roadmap2030\data\crop_maps\artichoke_harea.tif",
    #'asparagus': r"D:\repositories\roadmap2030\data\crop_maps\asparagus_harea.tif",
    #'avocado': r"D:\repositories\roadmap2030\data\crop_maps\avocado_harea.tif",
    #'bambara': r"D:\repositories\roadmap2030\data\crop_maps\bambara_harea.tif",
    #'banana': r"D:\repositories\roadmap2030\data\crop_maps\banana_harea.tif",
    #'barley': r"D:\repositories\roadmap2030\data\crop_maps\barley_harea.tif",
    #'bean': r"D:\repositories\roadmap2030\data\crop_maps\bean_harea.tif",
    #'beetfor': r"D:\repositories\roadmap2030\data\crop_maps\beetfor_harea.tif",
    #'berrynes': r"D:\repositories\roadmap2030\data\crop_maps\berrynes_harea.tif",
    #'blueberry': r"D:\repositories\roadmap2030\data\crop_maps\blueberry_harea.tif",
    #'brazil': r"D:\repositories\roadmap2030\data\crop_maps\brazil_harea.tif",
    #'broadbean': r"D:\repositories\roadmap2030\data\crop_maps\broadbean_harea.tif",
    #'buckwheat': r"D:\repositories\roadmap2030\data\crop_maps\buckwheat_harea.tif",
    #'cabbage': r"D:\repositories\roadmap2030\data\crop_maps\cabbage_harea.tif",
    #'cabbagefor': r"D:\repositories\roadmap2030\data\crop_maps\cabbagefor_harea.tif",
    #'canaryseed': r"D:\repositories\roadmap2030\data\crop_maps\canaryseed_harea.tif",
    #'carob': r"D:\repositories\roadmap2030\data\crop_maps\carob_harea.tif",
    #'carrot': r"D:\repositories\roadmap2030\data\crop_maps\carrot_harea.tif",
    #'carrotfor': r"D:\repositories\roadmap2030\data\crop_maps\carrotfor_harea.tif",
    #'cashew': r"D:\repositories\roadmap2030\data\crop_maps\cashew_harea.tif",
    #'cashewapple': r"D:\repositories\roadmap2030\data\crop_maps\cashewapple_harea.tif",
    #'cassava': r"D:\repositories\roadmap2030\data\crop_maps\cassava_harea.tif",
    #'castor': r"D:\repositories\roadmap2030\data\crop_maps\castor_harea.tif",
    #'cauliflower': r"D:\repositories\roadmap2030\data\crop_maps\cauliflower_harea.tif",
    #'cerealnes': r"D:\repositories\roadmap2030\data\crop_maps\cerealnes_harea.tif",
    #'cherry': r"D:\repositories\roadmap2030\data\crop_maps\cherry_harea.tif",
    #'chestnut': r"D:\repositories\roadmap2030\data\crop_maps\chestnut_harea.tif",
    #'chickpea': r"D:\repositories\roadmap2030\data\crop_maps\chickpea_harea.tif",
    #'chicory': r"D:\repositories\roadmap2030\data\crop_maps\chicory_harea.tif",
    #'chilleetc': r"D:\repositories\roadmap2030\data\crop_maps\chilleetc_harea.tif",
    #'cinnamon': r"D:\repositories\roadmap2030\data\crop_maps\cinnamon_harea.tif",
    #'citrusnes': r"D:\repositories\roadmap2030\data\crop_maps\citrusnes_harea.tif",
    #'clove': r"D:\repositories\roadmap2030\data\crop_maps\clove_harea.tif",
    #'clover': r"D:\repositories\roadmap2030\data\crop_maps\clover_harea.tif",
    #'cocoa': r"D:\repositories\roadmap2030\data\crop_maps\cocoa_harea.tif",
    #'coconut': r"D:\repositories\roadmap2030\data\crop_maps\coconut_harea.tif",
    #'coffee': r"D:\repositories\roadmap2030\data\crop_maps\coffee_harea.tif",
    #'coir': r"D:\repositories\roadmap2030\data\crop_maps\coir_harea.tif",
    #'cotton': r"D:\repositories\roadmap2030\data\crop_maps\cotton_harea.tif",
    #'cowpea': r"D:\repositories\roadmap2030\data\crop_maps\cowpea_harea.tif",
    #'cranberry': r"D:\repositories\roadmap2030\data\crop_maps\cranberry_harea.tif",
    #'cucumberetc': r"D:\repositories\roadmap2030\data\crop_maps\cucumberetc_harea.tif",
    #'currant': r"D:\repositories\roadmap2030\data\crop_maps\currant_harea.tif",
    #'date': r"D:\repositories\roadmap2030\data\crop_maps\date_harea.tif",
    #'eggplant': r"D:\repositories\roadmap2030\data\crop_maps\eggplant_harea.tif",
    #'fibrenes': r"D:\repositories\roadmap2030\data\crop_maps\fibrenes_harea.tif",
    #'fig': r"D:\repositories\roadmap2030\data\crop_maps\fig_harea.tif",
    #'flax': r"D:\repositories\roadmap2030\data\crop_maps\flax_harea.tif",
    #'fonio': r"D:\repositories\roadmap2030\data\crop_maps\fonio_harea.tif",
    #'fornes': r"D:\repositories\roadmap2030\data\crop_maps\fornes_harea.tif",
    #'fruitnes': r"D:\repositories\roadmap2030\data\crop_maps\fruitnes_harea.tif",
    #'garlic': r"D:\repositories\roadmap2030\data\crop_maps\garlic_harea.tif",
    #'ginger': r"D:\repositories\roadmap2030\data\crop_maps\ginger_harea.tif",
    #'gooseberry': r"D:\repositories\roadmap2030\data\crop_maps\gooseberry_harea.tif",
    #'grape': r"D:\repositories\roadmap2030\data\crop_maps\grape_harea.tif",
    #'grapefruitetc': r"D:\repositories\roadmap2030\data\crop_maps\grapefruitetc_harea.tif",
    #'grassnes': r"D:\repositories\roadmap2030\data\crop_maps\grassnes_harea.tif",
    #'greenbean': r"D:\repositories\roadmap2030\data\crop_maps\greenbean_harea.tif",
    #'greenbroadbean': r"D:\repositories\roadmap2030\data\crop_maps\greenbroadbean_harea.tif",
    #'greencorn': r"D:\repositories\roadmap2030\data\crop_maps\greencorn_harea.tif",
    #'greenonion': r"D:\repositories\roadmap2030\data\crop_maps\greenonion_harea.tif",
    #'greenpea': r"D:\repositories\roadmap2030\data\crop_maps\greenpea_harea.tif",
    #'groundnut': r"D:\repositories\roadmap2030\data\crop_maps\groundnut_harea.tif",
    #'gums': r"D:\repositories\roadmap2030\data\crop_maps\gums_harea.tif",
    #'hazelnut': r"D:\repositories\roadmap2030\data\crop_maps\hazelnut_harea.tif",
    #'hemp': r"D:\repositories\roadmap2030\data\crop_maps\hemp_harea.tif",
    #'hempseed': r"D:\repositories\roadmap2030\data\crop_maps\hempseed_harea.tif",
    #'hop': r"D:\repositories\roadmap2030\data\crop_maps\hop_harea.tif",
    #'jute': r"D:\repositories\roadmap2030\data\crop_maps\jute_harea.tif",
    #'jutelikefiber': r"D:\repositories\roadmap2030\data\crop_maps\jutelikefiber_harea.tif",
    #'kapokfiber': r"D:\repositories\roadmap2030\data\crop_maps\kapokfiber_harea.tif",
    #'kapokseed': r"D:\repositories\roadmap2030\data\crop_maps\kapokseed_harea.tif",
    #'karite': r"D:\repositories\roadmap2030\data\crop_maps\karite_harea.tif",
    #'kiwi': r"D:\repositories\roadmap2030\data\crop_maps\kiwi_harea.tif",
    #'kolanut': r"D:\repositories\roadmap2030\data\crop_maps\kolanut_harea.tif",
    #'legumenes': r"D:\repositories\roadmap2030\data\crop_maps\legumenes_harea.tif",
    #'lemonlime': r"D:\repositories\roadmap2030\data\crop_maps\lemonlime_harea.tif",
    #'lentil': r"D:\repositories\roadmap2030\data\crop_maps\lentil_harea.tif",
    #'lettuce': r"D:\repositories\roadmap2030\data\crop_maps\lettuce_harea.tif",
    #'linseed': r"D:\repositories\roadmap2030\data\crop_maps\linseed_harea.tif",
    #'lupin': r"D:\repositories\roadmap2030\data\crop_maps\lupin_harea.tif",
    #'maize': r"D:\repositories\roadmap2030\data\crop_maps\maize_harea.tif",
    #'maizefor': r"D:\repositories\roadmap2030\data\crop_maps\maizefor_harea.tif",
    #'mango': r"D:\repositories\roadmap2030\data\crop_maps\mango_harea.tif",
    #'mate': r"D:\repositories\roadmap2030\data\crop_maps\mate_harea.tif",
    #'melonetc': r"D:\repositories\roadmap2030\data\crop_maps\melonetc_harea.tif",
    #'melonseed': r"D:\repositories\roadmap2030\data\crop_maps\melonseed_harea.tif",
    #'millet': r"D:\repositories\roadmap2030\data\crop_maps\millet_harea.tif",
    #'mixedgrain': r"D:\repositories\roadmap2030\data\crop_maps\mixedgrain_harea.tif",
    #'mixedgrass': r"D:\repositories\roadmap2030\data\crop_maps\mixedgrass_harea.tif",
    #'mushroom': r"D:\repositories\roadmap2030\data\crop_maps\mushroom_harea.tif",
    #'mustard': r"D:\repositories\roadmap2030\data\crop_maps\mustard_harea.tif",
    #'nutmeg': r"D:\repositories\roadmap2030\data\crop_maps\nutmeg_harea.tif",
    #'nutnes': r"D:\repositories\roadmap2030\data\crop_maps\nutnes_harea.tif",
    #'oats': r"D:\repositories\roadmap2030\data\crop_maps\oats_harea.tif",
    #'oilpalm': r"D:\repositories\roadmap2030\data\crop_maps\oilpalm_harea.tif",
    #'oilseedfor': r"D:\repositories\roadmap2030\data\crop_maps\oilseedfor_harea.tif",
    #'oilseednes': r"D:\repositories\roadmap2030\data\crop_maps\oilseednes_harea.tif",
    #'okra': r"D:\repositories\roadmap2030\data\crop_maps\okra_harea.tif",
    #'olive': r"D:\repositories\roadmap2030\data\crop_maps\olive_harea.tif",
    #'onion': r"D:\repositories\roadmap2030\data\crop_maps\onion_harea.tif",
    #'orange': r"D:\repositories\roadmap2030\data\crop_maps\orange_harea.tif",
    #'papaya': r"D:\repositories\roadmap2030\data\crop_maps\papaya_harea.tif",
    #'pea': r"D:\repositories\roadmap2030\data\crop_maps\pea_harea.tif",
    #'peachetc': r"D:\repositories\roadmap2030\data\crop_maps\peachetc_harea.tif",
    #'pear': r"D:\repositories\roadmap2030\data\crop_maps\pear_harea.tif",
    #'pepper': r"D:\repositories\roadmap2030\data\crop_maps\pepper_harea.tif",
    #'peppermint': r"D:\repositories\roadmap2030\data\crop_maps\peppermint_harea.tif",
    #'persimmon': r"D:\repositories\roadmap2030\data\crop_maps\persimmon_harea.tif",
    #'pigeonpea': r"D:\repositories\roadmap2030\data\crop_maps\pigeonpea_harea.tif",
    #'pimento': r"D:\repositories\roadmap2030\data\crop_maps\pimento_harea.tif",
    #'pineapple': r"D:\repositories\roadmap2030\data\crop_maps\pineapple_harea.tif",
    #'pistachio': r"D:\repositories\roadmap2030\data\crop_maps\pistachio_harea.tif",
    #'plantain': r"D:\repositories\roadmap2030\data\crop_maps\plantain_harea.tif",
    #'plum': r"D:\repositories\roadmap2030\data\crop_maps\plum_harea.tif",
    #'popcorn': r"D:\repositories\roadmap2030\data\crop_maps\popcorn_harea.tif",
    #'poppy': r"D:\repositories\roadmap2030\data\crop_maps\poppy_harea.tif",
    #'potato': r"D:\repositories\roadmap2030\data\crop_maps\potato_harea.tif",
    #'pulsenes': r"D:\repositories\roadmap2030\data\crop_maps\pulsenes_harea.tif",
    #'pumpkinetc': r"D:\repositories\roadmap2030\data\crop_maps\pumpkinetc_harea.tif",
    #'pyrethrum': r"D:\repositories\roadmap2030\data\crop_maps\pyrethrum_harea.tif",
    #'quince': r"D:\repositories\roadmap2030\data\crop_maps\quince_harea.tif",
    #'quinoa': r"D:\repositories\roadmap2030\data\crop_maps\quinoa_harea.tif",
    #'ramie': r"D:\repositories\roadmap2030\data\crop_maps\ramie_harea.tif",
    #'rapeseed': r"D:\repositories\roadmap2030\data\crop_maps\rapeseed_harea.tif",
    #'rasberry': r"D:\repositories\roadmap2030\data\crop_maps\rasberry_harea.tif",
    #'rice': r"D:\repositories\roadmap2030\data\crop_maps\rice_harea.tif",
    #'rootnes': r"D:\repositories\roadmap2030\data\crop_maps\rootnes_harea.tif",
    #'rubber': r"D:\repositories\roadmap2030\data\crop_maps\rubber_harea.tif",
    #'rye': r"D:\repositories\roadmap2030\data\crop_maps\rye_harea.tif",
    #'ryefor': r"D:\repositories\roadmap2030\data\crop_maps\ryefor_harea.tif",
    #'safflower': r"D:\repositories\roadmap2030\data\crop_maps\safflower_harea.tif",
    #'sesame': r"D:\repositories\roadmap2030\data\crop_maps\sesame_harea.tif",
    #'sisal': r"D:\repositories\roadmap2030\data\crop_maps\sisal_harea.tif",
    #'sorghum': r"D:\repositories\roadmap2030\data\crop_maps\sorghum_harea.tif",
    #'sorghumfor': r"D:\repositories\roadmap2030\data\crop_maps\sorghumfor_harea.tif",
    #'sourcherry': r"D:\repositories\roadmap2030\data\crop_maps\sourcherry_harea.tif",
    #'soybean': r"D:\repositories\roadmap2030\data\crop_maps\soybean_harea.tif",
    #'spicenes': r"D:\repositories\roadmap2030\data\crop_maps\spicenes_harea.tif",
    #'spinach': r"D:\repositories\roadmap2030\data\crop_maps\spinach_harea.tif",
    #'stonefruitnes': r"D:\repositories\roadmap2030\data\crop_maps\stonefruitnes_harea.tif",
    #'strawberry': r"D:\repositories\roadmap2030\data\crop_maps\strawberry_harea.tif",
    #'stringbean': r"D:\repositories\roadmap2030\data\crop_maps\stringbean_harea.tif",
    #'sugarbeet': r"D:\repositories\roadmap2030\data\crop_maps\sugarbeet_harea.tif",
    #'sugarcane': r"D:\repositories\roadmap2030\data\crop_maps\sugarcane_harea.tif",
    #'sugarnes': r"D:\repositories\roadmap2030\data\crop_maps\sugarnes_harea.tif",
    #'sunflower': r"D:\repositories\roadmap2030\data\crop_maps\sunflower_harea.tif",
    #'swedefor': r"D:\repositories\roadmap2030\data\crop_maps\swedefor_harea.tif",
    #'sweetpotato': r"D:\repositories\roadmap2030\data\crop_maps\sweetpotato_harea.tif",
    #'tangetc': r"D:\repositories\roadmap2030\data\crop_maps\tangetc_harea.tif",
    #'taro': r"D:\repositories\roadmap2030\data\crop_maps\taro_harea.tif",
    #'tea': r"D:\repositories\roadmap2030\data\crop_maps\tea_harea.tif",
    #'tobacco': r"D:\repositories\roadmap2030\data\crop_maps\tobacco_harea.tif",
    #'tomato': r"D:\repositories\roadmap2030\data\crop_maps\tomato_harea.tif",
    #'triticale': r"D:\repositories\roadmap2030\data\crop_maps\triticale_harea.tif",
    #'tropicalnes': r"D:\repositories\roadmap2030\data\crop_maps\tropicalnes_harea.tif",
    #'tung': r"D:\repositories\roadmap2030\data\crop_maps\tung_harea.tif",
    #'turnipfor': r"D:\repositories\roadmap2030\data\crop_maps\turnipfor_harea.tif",
    #'vanilla': r"D:\repositories\roadmap2030\data\crop_maps\vanilla_harea.tif",
    #'vegetablenes': r"D:\repositories\roadmap2030\data\crop_maps\vegetablenes_harea.tif",
    #'vegfor': r"D:\repositories\roadmap2030\data\crop_maps\vegfor_harea.tif",
    #'vetch': r"D:\repositories\roadmap2030\data\crop_maps\vetch_harea.tif",
    #'walnut': r"D:\repositories\roadmap2030\data\crop_maps\walnut_harea.tif",
    #'watermelon': r"D:\repositories\roadmap2030\data\crop_maps\watermelon_harea.tif",
    #'wheat': r"D:\repositories\roadmap2030\data\crop_maps\wheat_harea.tif",
    #'yam': r"D:\repositories\roadmap2030\data\crop_maps\yam_harea.tif",
    #'yautia': r"D:\repositories\roadmap2030\data\crop_maps\yautia_harea.tif",
    #'2020_QF': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_2020_results\amazon_swy_2020_results_QF.tif",
    #'2020_QF_1': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_2020_results\amazon_swy_2020_results_QF_1.tif",
    #'2020_QF_2': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_2020_results\amazon_swy_2020_results_QF_2.tif",
    #'2020_QF_3': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_2020_results\amazon_swy_2020_results_QF_3.tif",
    #'2020_QF_4': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_2020_results\amazon_swy_2020_results_QF_4.tif",
    #'2020_QF_5': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_2020_results\amazon_swy_2020_results_QF_5.tif",
    #'2020_QF_6': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_2020_results\amazon_swy_2020_results_QF_6.tif",
    #'2020_QF_7': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_2020_results\amazon_swy_2020_results_QF_7.tif",
    #'2020_QF_8': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_2020_results\amazon_swy_2020_results_QF_8.tif",
    #'2020_QF_9': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_2020_results\amazon_swy_2020_results_QF_9.tif",
    #'2020_QF_10': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_2020_results\amazon_swy_2020_results_QF_10.tif",
    #'2020_QF_11': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_2020_results\amazon_swy_2020_results_QF_11.tif",
    #'2020_QF_12': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_2020_results\amazon_swy_2020_results_QF_12.tif",
    #'2020_B_sum': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_2020_results\amazon_swy_2020_results_B_sum.tif",
    #'1992_QF': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_1992_results\QF.tif",
    #'1992_QF_1': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_1992_results\QF_1.tif",
    #'1992_QF_2': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_1992_results\QF_2.tif",
    #'1992_QF_3': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_1992_results\QF_3.tif",
    #'1992_QF_4': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_1992_results\QF_4.tif",
    #'1992_QF_5': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_1992_results\QF_5.tif",
    #'1992_QF_6': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_1992_results\QF_6.tif",
    #'1992_QF_7': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_1992_results\QF_7.tif",
    #'1992_QF_8': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_1992_results\QF_8.tif",
    #'1992_QF_9': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_1992_results\QF_9.tif",
    #'1992_QF_10': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_1992_results\QF_10.tif",
    #'1992_QF_11': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_1992_results\QF_11.tif",
    #'1992_QF_12': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_1992_results\QF_12.tif",
    #'1992_B_sum': r"D:\repositories\roadmap2030\global_swy_runs\amazon_swy_1992_results\B_sum.tif",
    #'ARPA_AGB_from_Karuna_CO2e': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\ARPA_AGB_from_Karuna_CO2e_compressed.tif",
    #'carbon_2008': r"D:\repositories\roadmap2030\data\ABUNCHASERVICES\fc_stack_hansen_forest_cover2008_compressed_std_forest_edge_result.tif",
}


VECTOR_PATH_LOOKUP = {
    #'37_GEF_Peru': r"D:\repositories\roadmap2030\data\37_GEF_Peru.gpkg",
    #'49_GEF_Guyana_KPMA_NRW': r"D:\repositories\roadmap2030\data\49_GEF_Guyana_KPMA_NRW.gpkg",
    #'105_Arpa_nonoverlapping_clipped': r"D:\repositories\roadmap2030\data\105_Arpa_nonoverlapping_clipped.gpkg",
    #'106_HECO': r"D:\repositories\roadmap2030\data\106_HECO.gpkg",
    #'110_PatriomonioPeru': r"D:\repositories\roadmap2030\data\110_PatriomonioPeru.gpkg",
    #'118_NBSOP2': r"D:\repositories\roadmap2030\data\118_NBSOP2.gpkg",
    #'118_NBSOP3': r"D:\repositories\roadmap2030\data\118_NBSOP3.gpkg",
    #'118_NBSOP4': r"D:\repositories\roadmap2030\data\118_NBSOP4.gpkg",
    #'120_Tapajos': r"D:\repositories\roadmap2030\data\120_Tapajos.gpkg",
    #'133_Sall': r"D:\repositories\roadmap2030\data\133_Sall.gpkg",
    #'All_WWF_Amazon': r"D:\repositories\roadmap2030\data\000_All_WWF_Amazon_sites.gpkg",
    #'Amazon': r"D:\repositories\roadmap2030\data\Amazon_aoi.gpkg",
    #'non-arpa': r"D:\repositories\roadmap2030\data\non-Arpa_nonoverlapping-in-m.gpkg",
    #'arpa': r"D:\repositories\roadmap2030\data\Arpa_nonoverlapping-in-m.gpkg",
    #'colombia': r"D:\repositories\roadmap2030\data\Colombia.gpkg",
    #'peru': r"D:\repositories\roadmap2030\data\Peru.gpkg",
    #'tapajos': r"D:\repositories\roadmap2030\data\Tapajos.gpkg",
    #'NGP': r"D:\repositories\roadmap2030\data\NGP.gpkg",
    #'RGBR': r"D:\repositories\roadmap2030\data\RGRB.gpkg",
    #'Arctic': r"D:\repositories\roadmap2030\data\Arctic.gpkg",
    #'105': r"D:\repositories\roadmap2030\data\aois\final_pilot\105.gpkg",
    #'106': r"D:\repositories\roadmap2030\data\aois\final_pilot\106.gpkg",
    #'110': r"D:\repositories\roadmap2030\data\aois\final_pilot\110.gpkg",
    #'118_combined': r"D:\repositories\roadmap2030\data\aois\final_pilot\118_combined.gpkg",
    #'120': r"D:\repositories\roadmap2030\data\aois\final_pilot\120.gpkg",
    #'133': r"D:\repositories\roadmap2030\data\aois\final_pilot\133.gpkg",
    #'199': r"D:\repositories\roadmap2030\data\aois\final_pilot\158-199.gpkg",
    #'196': r"D:\repositories\roadmap2030\data\aois\final_pilot\196.gpkg",
    #'23': r"D:\repositories\roadmap2030\data\aois\final_pilot\23.gpkg",
    #'312': r"D:\repositories\roadmap2030\data\aois\final_pilot\292-299-312.gpkg",
    #'312_wg84': r"D:\repositories\roadmap2030\data\aois\final_pilot\292-299-312_wg84.gpkg",
    #'302': r"D:\repositories\roadmap2030\data\aois\final_pilot\300-302.gpkg",
    #'313': r"D:\repositories\roadmap2030\data\aois\final_pilot\313.gpkg",
    #'37': r"D:\repositories\roadmap2030\data\aois\final_pilot\37.gpkg",
    #'49': r"D:\repositories\roadmap2030\data\aois\final_pilot\49.gpkg",
    'costa_rica_pfp': r"D:\repositories\roadmap2030\data\aois\costa_rica_pfp_in_m.gpkg",
    'costa_rica': r"D:\repositories\roadmap2030\data\aois\costa_rica_in_m.gpkg",
    'great_bear_pfp': r"D:\repositories\roadmap2030\data\aois\great_bear_pfp_in_m.gpkg",
    'british_columbia': r"D:\repositories\roadmap2030\data\aois\british_columbia_in_m.gpkg"
}

OUTPUT_DIR = './results'
CLIPPED_DIR = os.path.join(OUTPUT_DIR, 'clipped')
for dirpath in [OUTPUT_DIR, CLIPPED_DIR]:
    os.makedirs(dirpath, exist_ok=True)


def vector_area_in_ha(vector_path):
    dataset = gdal.OpenEx(vector_path, gdal.OF_VECTOR)
    layer = dataset.GetLayer()

    source_srs = layer.GetSpatialRef()
    target_srs = osr.SpatialReference()
    target_srs.ImportFromProj4("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +R=6371007 +units=m +no_defs")

    total_area_m2 = 0

    for feature in layer:
        geometry = feature.GetGeometryRef()
        if geometry is not None:
            geom_clone = geometry.Clone()
            geom_clone.AssignSpatialReference(source_srs)
            geom_clone.TransformTo(target_srs)
            total_area_m2 += geom_clone.GetArea()

    total_area_ha = total_area_m2 / 10000
    return total_area_ha


def create_subset(gdf, name, target_vector_path):
    LOGGER.info(f'creating subset of {name}')
    subset_gdf = gdf[gdf["Name"] == name]
    subset_gdf.to_file(target_vector_path, driver="GPKG")
    LOGGER.info(f'done with subset of {name}')


def clip_raster(base_raster_path, summary_vector_path, temp_clip_path):
    base_raster_info = geoprocessing.get_raster_info(base_raster_path)
    summary_vector_info = geoprocessing.get_vector_info(summary_vector_path)
    target_pixel_size = base_raster_info['pixel_size']
    base_vector_bb = summary_vector_info['bounding_box']

    target_bb = geoprocessing.transform_bounding_box(
        base_vector_bb, summary_vector_info['projection_wkt'],
        base_raster_info['projection_wkt'])

    geoprocessing.warp_raster(
        base_raster_path, target_pixel_size, temp_clip_path,
        'near', target_bb=target_bb, vector_mask_options={
            'mask_vector_path': summary_vector_path,
            'all_touched': True})


def get_area_stats(raster_path, thresholds):
    r = gdal.Open(raster_path)
    proj_wkt = r.GetProjection()
    srs = osr.SpatialReference()
    srs.ImportFromWkt(proj_wkt)
    if srs.IsGeographic():
        warp_srs = osr.SpatialReference()
        warp_srs.ImportFromEPSG(54009)  # Mollweide (World)
        warped_ds = gdal.Warp(
            '', r, format='MEM',
            dstSRS=warp_srs,
            resampleAlg=gdal.GRA_NearestNeighbour
        )
    else:
        warped_ds = r

    band = warped_ds.GetRasterBand(1)
    arr = band.ReadAsArray().astype(float)
    nodata = band.GetNoDataValue()

    valid_mask = numpy.ones_like(arr, dtype=bool)
    if nodata is not None:
        valid_mask &= (arr != nodata)
    valid_mask &= ~numpy.isnan(arr)
    arr = arr[valid_mask]

    gt = warped_ds.GetGeoTransform()
    pixel_area_m2 = abs(gt[1] * gt[5])
    results = {}

    for thr in thresholds:
        pix_count = numpy.count_nonzero(arr >= thr)
        area_ha = (pix_count * pixel_area_m2) / 10000.0  # 1 ha = 10,000 m^2
        results[f'area_ge_{thr}'] = area_ha

    area_ha = (arr.size * pixel_area_m2) / 10000.0  # 1 ha = 10,000 m^2
    results['area_ha'] = area_ha

    return results


def get_stats(raster_path):
    r = gdal.OpenEx(raster_path)
    b = r.GetRasterBand(1)
    array = b.ReadAsArray()
    nodata = b.GetNoDataValue()
    array = array[(array != nodata) & ~numpy.isnan(array)]

    stats = {
        'min': numpy.min(array),
        'max': numpy.max(array),
        'sum': numpy.sum(array),
        'mean': numpy.mean(array)
    }

    percentile_dict = {
        f'p{percentile}': value
        for percentile, value in zip(
            PERCENTILES_LIST,
            numpy.percentile(array, PERCENTILES_LIST))
    }
    stats.update(percentile_dict)

    value_thresholds = get_area_stats(raster_path, THRESHOLD_AREA_LIST)

    stats.update(value_thresholds)


    return stats


def dump_results_to_csv(results, vector_path_lookup, csv_path):
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(
            ["vector_id", "raster_name", "area_ha", "min", "max", "mean", "sum"] +
            [f'p{percentile}' for percentile in PERCENTILES_LIST] +
            [f'area_ge_{threshold}' for threshold in THRESHOLD_AREA_LIST])
        for vector_id, info_dict in results.items():
            area_ha = info_dict.get("area_ha", None)
            if area_ha is not None:
                writer.writerow([
                    vector_id,
                    "",         # raster_name is empty
                    area_ha,
                    "", "", "", ""  # no min/max/mean/sum for area
                ])

            # Next rows: stats for each raster
            # (anything that's not "area_ha" in `results[vector_id]`)
            for raster_basename, stats_dict in info_dict.items():
                if raster_basename == "area_ha":
                    continue
                if not isinstance(stats_dict, dict):
                    continue

                # Extract stats if they exist
                r_min = stats_dict.get("min", "")
                r_max = stats_dict.get("max", "")
                r_mean = stats_dict.get("mean", "")
                r_sum = stats_dict.get("sum", "")
                area_ha = stats_dict.get('area_ha', '')
                writer.writerow([
                    vector_id,
                    raster_basename,  # raster_name
                    area_ha,               # area_ha is empty here
                    r_min,
                    r_max,
                    r_mean,
                    r_sum
                ] + [stats_dict.get(f'p{percentile}', '') for percentile in PERCENTILES_LIST] +
                [stats_dict.get(f'area_ge_{thr}', '') for thr in THRESHOLD_AREA_LIST])


def main():
    """Entry point."""
    print(os.cpu_count())
    task_graph = taskgraph.TaskGraph(OUTPUT_DIR, os.cpu_count(), reporting_interval=10.0)
    results = collections.defaultdict(lambda: collections.defaultdict(dict))
    for vector_id, vector_path in VECTOR_PATH_LOOKUP.items():
        results[vector_id]['area_ha'] = vector_area_in_ha(vector_path)
        LOGGER.info(f'processing {vector_id}')
        for raster_basename, raster_path in BASE_RASTER_LOOKUP.items():
            if not os.path.exists(raster_path):
                raise RuntimeError(f'{raster_path} not found')
            LOGGER.info(f'clipping {raster_basename} to {vector_id}')
            clipped_raster_path = os.path.join(
                CLIPPED_DIR, f'{vector_id}_{raster_basename}.tif')
            clipped_task = task_graph.add_task(
                func=clip_raster,
                args=(raster_path, vector_path, clipped_raster_path),
                ignore_path_list=[vector_path],
                target_path_list=[clipped_raster_path],
                task_name=f'clipping {raster_path} to {vector_path}')

            stats_task = task_graph.add_task(
                func=get_stats,
                args=(clipped_raster_path,),
                dependent_task_list=[clipped_task],
                store_result=True,
                task_name=f'stats for {raster_path}')
            results[vector_id][raster_basename]['stats'] = stats_task

    for vector_id in VECTOR_PATH_LOOKUP:
        for raster_basename in BASE_RASTER_LOOKUP:
            stats_task = results[vector_id][raster_basename]['stats']
            results[vector_id][raster_basename] = {}
            for fn_str, value in stats_task.get().items():
                results[vector_id][raster_basename][fn_str] = value

    task_graph.join()
    print(results)
    timestamp = datetime.datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    output_filename = f'results_{timestamp}.csv'
    dump_results_to_csv(results, VECTOR_PATH_LOOKUP, output_filename)
    print(f'all done -- results in {output_filename}!')


if __name__ == '__main__':
    main()
