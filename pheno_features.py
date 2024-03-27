# pheno_features_definitions.py
#######################################
############ PFT features #############
#######################################

# Correct lung function features for confounders
PFT_FEATURES = ['Resting_SaO2_P2' ] ##features
PFT_NAMES = [r'$\mathbf{Resting\ SaO_{2}}$'] ## names to be printed on plot
PFT_order_names = ['Age_P2', 'C(gender)[T.2]','C(race)[T.2]','Height_CM_P2',
            "C(Disease, Treatment(reference='Control'))[T.Asthma]",
            "C(Disease, Treatment(reference='Control'))[T.COPD]", 
            "C(Disease, Treatment(reference='Control'))[T.ACO]"
            ]
PFT_order_names = PFT_order_names[::-1]

PFT_rename = {'Age_P2':'Age', 
        'C(gender)[T.2]':'Gender',
        'C(race)[T.2]':'Race',
        'Height_CM_P2':'Height',  
        "C(Disease, Treatment(reference='Control'))[T.Asthma]": 'Asthma', 
        "C(Disease, Treatment(reference='Control'))[T.COPD]":'COPD', 
        "C(Disease, Treatment(reference='Control'))[T.ACO]":'ACO'}

#######################################
############ CT features ##############
#######################################

            
# Define outcomes
CT_FEATURES = [
    'pctEmph_Thirona_P2',
    'pctGasTrap_Thirona_P2',
    'Pi10_Thirona_P2',
    'Perc15_density_Thirona_P2', 
    "PRM_pct_airtrapping_Thirona_P2",
    'WallAreaPct_seg_Thirona_P2',
  ]

# Define names corresponding to each feature
CT_NAMES = [
    '% Emphysema (<-950HU)','Gas Trapping (<-856HU)', 'Pi10','Lung density',
    r'$\mathbf{PRM}^{fSAD}$',
    'Wall Area (%)',  
 
]

CT_order_names = [
      'Age_P2',
       'BMI_P2','C(gender)[T.2]','C(race)[T.2]', 'C(smoking_status_P2)[T.2.0]',  # 'gender', 'race',
     "C(Disease, Treatment(reference='Control'))[T.Asthma]",
       "C(Disease, Treatment(reference='Control'))[T.COPD]",
     "C(Disease, Treatment(reference='Control'))[T.ACO]",
    'C(scannerId_P2)[T.C0201]', 'C(scannerId_P2)[T.C0301]',
       'C(scannerId_P2)[T.C0302]', 'C(scannerId_P2)[T.C0402]',
       'C(scannerId_P2)[T.C0403]', 'C(scannerId_P2)[T.C0501]',
       'C(scannerId_P2)[T.C0502]', 'C(scannerId_P2)[T.C0503]',
       'C(scannerId_P2)[T.C0602]', 'C(scannerId_P2)[T.C0603]',
       'C(scannerId_P2)[T.C0701]', 'C(scannerId_P2)[T.C0703]',
       'C(scannerId_P2)[T.C0806]', 'C(scannerId_P2)[T.C0807]',
       'C(scannerId_P2)[T.C0902]', 'C(scannerId_P2)[T.C0904]',
       'C(scannerId_P2)[T.C1004]', 'C(scannerId_P2)[T.C1102]',
       'C(scannerId_P2)[T.C1103]', 'C(scannerId_P2)[T.C1202]',
       'C(scannerId_P2)[T.C1302]', 'C(scannerId_P2)[T.C1402]',
       'C(scannerId_P2)[T.C1502]', 'C(scannerId_P2)[T.C1503]',
       'C(scannerId_P2)[T.C1605]', 'C(scannerId_P2)[T.C1704]',
       'C(scannerId_P2)[T.C1707]', 'C(scannerId_P2)[T.C1708]',
       'C(scannerId_P2)[T.C1901]', 'C(scannerId_P2)[T.C1903]',
       'C(scannerId_P2)[T.C2001]', 'C(scannerId_P2)[T.C2103]'
]
CT_order_names = CT_order_names[::-1]

CT_rename = {'Age_P2':'Age', 
          'C(gender)[T.2]':'Gender',
          'C(race)[T.2]':'Race',
          'Height_CM_P2':'Height',  
          'BMI_P2': 'BMI',
          'C(smoking_status_P2)[T.2.0]':'Smoking Status',
          'smoking_status_P2':'Smoking Status' ,
          "C(Disease, Treatment(reference='Control'))[T.Asthma]": 'Asthma', 
          "C(Disease, Treatment(reference='Control'))[T.COPD]":'COPD', 
           "C(Disease, Treatment(reference='Control'))[T.ACO]":'ACO',
          'C(scannerId_P2)[T.C0201]':'C0201',
          'C(scannerId_P2)[T.C0301]':'C0301',
       'C(scannerId_P2)[T.C0302]': 'C0302',
        'C(scannerId_P2)[T.C0402]':'C0402',
       'C(scannerId_P2)[T.C0403]':'C0403', 
        'C(scannerId_P2)[T.C0501]':'C0501',
       'C(scannerId_P2)[T.C0502]':'C0502',
       'C(scannerId_P2)[T.C0503]':'C0503',
       'C(scannerId_P2)[T.C0602]':'C0602', 
          'C(scannerId_P2)[T.C0603]':'C0603',
       'C(scannerId_P2)[T.C0701]':'C0701', 
          'C(scannerId_P2)[T.C0703]':'C0703',
       'C(scannerId_P2)[T.C0806]':'C0806', 
          'C(scannerId_P2)[T.C0807]':'C0807',
       'C(scannerId_P2)[T.C0902]':'C0902', 
          'C(scannerId_P2)[T.C0904]':'C0904',
       'C(scannerId_P2)[T.C1004]':'C1004', 
          'C(scannerId_P2)[T.C1102]':'C1102',
       'C(scannerId_P2)[T.C1103]':'C1103', 
          'C(scannerId_P2)[T.C1202]':'C1202',
       'C(scannerId_P2)[T.C1302]':'C1302', 
          'C(scannerId_P2)[T.C1402]':'C1402',
       'C(scannerId_P2)[T.C1502]':'C1502', 
          'C(scannerId_P2)[T.C1503]':'C1503',
       'C(scannerId_P2)[T.C1605]':'C1605', 
          'C(scannerId_P2)[T.C1704]':'C1704',
       'C(scannerId_P2)[T.C1707]':'C1707', 
          'C(scannerId_P2)[T.C1708]':'C1708',
       'C(scannerId_P2)[T.C1901]':'C1901', 
          'C(scannerId_P2)[T.C1903]':'C1903',
       'C(scannerId_P2)[T.C2001]':'C2001', 
          'C(scannerId_P2)[T.C2103]':'C2103'
         }



