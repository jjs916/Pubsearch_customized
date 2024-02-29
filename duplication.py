import pandas as pd
from collections import Counter

# Bring File
file_path = 'DBLPpub_post_processed.xlsx'
sheet_name = 'database'
data = pd.read_excel(file_path, sheet_name=sheet_name)

# Set year range and target publication titles
year_range = range(2019, 2024)
target_journals = [#IMPORTANT: Every Publication titles that needed to be check should be in this list!
    'ICRA', 'IROS', 'IEEE Robotics Autom. Lett.', 'IEEE Trans. Robotics',
    'Int. J. Robotics Res.', 'Sci. Robotics', 'Robotics: Science and Systems',
    'CoRL', 'RO-MAN', 'CDC', 'ACC', 'ECC', 'IEEE Control. Syst. Lett.',
    'IEEE Trans. Autom. Control.', 'CVPR', 'NeurIPS'
]

university_list = [#IMPORTANT: Every Institutions that needed to be check should be in this list!
    'BHT Berlin', 'Constructor Uni Bremen', 'DLR', 'FAU Erlangen', 'Fraunhofer',
    'FU Berlin', 'HU Berlin', 'KIT', 'LMU', 'LU Hannover', 'MPI', 'RWTH Aachen', 'TU Berlin',
    'TU Braunschweig', 'TU Chemnitz', 'TU Darmstadt', 'TU Dortmund', 'TU Dresden', 'TU Ilmenau',
    'TU Kaiserslautern', 'TUM', 'Uni Augsburg', 'Uni Bamberg', 'Uni Bayreuth', 'Uni Bielefeld',
    'Uni Bonn', 'Uni Bremen', 'Uni Freiburg', 'Uni Hamburg', 'Uni Heidelberg', 'Uni Magdeburg',
    'Uni Stuttgart', 'Uni Tübingen', 'Universität der Künste Berlin', 'Universitätsklinikums Hamburg-Eppendorf', 'University of Oldenburg',
    'University of Osnabrück', 'UTN'
]


# Select second and third column
selected_rows = data[
    (data.iloc[:, 1].isin(year_range)) &
    (data.iloc[:, 2].isin(target_journals))
]

# Find duplicated titles
duplicated_indices = selected_rows[selected_rows.duplicated(subset=[selected_rows.columns[3], selected_rows.columns[4]], keep=False)].index
duplicated_rows = selected_rows.loc[duplicated_indices].sort_values(by=selected_rows.columns[3])


# print(selected_rows.duplicated(subset=selected_rows.columns[3], keep=False))

grouped = duplicated_rows.groupby(duplicated_rows.columns[3])

# Count Affiliation
result = {}
for name, group in grouped:
    university_counts = Counter(group.iloc[:, 4])
    result[name] = university_counts

university_journal_table = pd.DataFrame(0, index=university_list, columns=target_journals)

#Count duplicated numbers
for name, frequencies in result.items():
    for university, count in frequencies.items():
        group_with_university = grouped.get_group(name)
        third_column_values = group_with_university.iloc[:, 2].unique()
        if count >= 2 and university in university_list:
            for third_value in third_column_values:
                university_journal_table.loc[university, third_value] += (count-1)

output_file_path = 'duplication.xlsx'
university_journal_table.to_excel(output_file_path)
print(f"Data has been saved to {output_file_path}")