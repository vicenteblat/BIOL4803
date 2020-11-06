import argparse
import pandas

parser = argparse.ArgumentParser(description = "This program analyzes data taken from competition experiments to compare fitness between two strains")
parser.add_argument("ExcelFile", help = "Excel file contataining the competition experiments data")
parser.add_argument("OutputFile", help = "String representation of the desired name for the output file")
args = parser.parse_args()

well_data = pandas.read_excel(args.ExcelFile, sheet_name = 'WellData')
well_data = well_data.sort_values(by=['Well', 'Concentration'])

replicates = pandas.read_excel(args.ExcelFile, sheet_name = 'Dictionary')
replicates = replicates.sort_values(by=['Replicate','Time Point'])

output_df_1 = pandas.DataFrame(columns = ['Strain1','Strain2','Replicate','2','4','6'])

for i, row in replicates.iterrows():
    well = row['Well']
    time_point = row['Time Point']

    if time_point == 2:
        strain1 = row['Strain1']
        strain2 = row['Strain2']
        replicate = row['Replicate']
        low = well_data[well_data.Well == well].Concentration.values[0]
        high = well_data[well_data.Well == well].Concentration.values[1]
        value = low / (low + high)
        temp_row = pandas.DataFrame({'Strain1':strain1, 'Strain2':strain2, 'Replicate':replicate, '2':[value], '4':[0], '6':[0]})
    elif time_point == 4:
        low = well_data[well_data.Well == well].Concentration.values[0]
        high = well_data[well_data.Well == well].Concentration.values[1]
        value = low / (low + high)
        temp_row['4'] = value
    else:
        low = well_data[well_data.Well == well].Concentration.values[0]
        high = well_data[well_data.Well == well].Concentration.values[1]
        value = low / (low + high)
        temp_row['6'] = value
        output_df_1 = output_df_1.append(temp_row)

output_df_2 = pandas.DataFrame({'Strain1':strain1, 'Strain2':strain2, '2':[output_df_1['2'].mean()], '4':[output_df_1['4'].mean()], '6':[output_df_1['6'].mean()]})

writer = pandas.ExcelWriter(args.OutputFile + '.xlsx')
output_df_1.to_excel(writer, 'All_data', index = False, merge_cells=True)
output_df_2.to_excel(writer, 'Average_data', index = False, merge_cells=True)
writer.save()
