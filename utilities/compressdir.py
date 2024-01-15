import pandas as pd
import os

# replace with your folder's path
folder_path = 'pot_output\qvalue\output' #'C:\Users\ashworjs\OneDrive - Rose-Hulman Institute of Technology\Desktop\OOPS\pot_output\hq1\output'

all_files = os.listdir(folder_path)

# Filter out non-CSV files
csv_files = [f for f in all_files if f.endswith('.csv')]

print('a')
print(csv_files)

# Create a list to hold the dataframes
df_list = []

for csv in csv_files:
    file_path = os.path.join(folder_path, csv)
    df = pd.read_csv(file_path)
    df_list.append(df)

# Concatenate all data into one DataFrame
big_df = pd.concat(df_list, ignore_index=True)

# Save the final result to a new CSV file
big_df.to_csv(os.path.join(folder_path, 'Qvalue.csv'), index=False)
