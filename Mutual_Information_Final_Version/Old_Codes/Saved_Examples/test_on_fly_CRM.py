import os
import Mutual_Information_Final_Version.functions.mutual_information_functions.cc_functions as pcmi
os.chdir("../../../")
raw_data_folder_name = f"foxo_igfr_3"
raw_data_folder_location = "Mutual_Information_Final_Version/Data/IGFR/Raw_Cell_Files_Conditioning/"
pcmi_value, cc_in = pcmi.pcmi_at_cc_generate_CRM_on_fly(raw_data_folder_name, dir=raw_data_folder_location,lock_parameter={"time":"90_min"},assumed_distribution="gamma")