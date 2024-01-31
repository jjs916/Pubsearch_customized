from scholarly import scholarly
import csv
import pandas as pd
import re

# Read Excel file and extract author names
file_path = './PIlist_MIRMI.xlsx'  # Replace with your file path
df = pd.read_excel(file_path)
# Assuming the column containing author names is named 'Authors'
PI = df['PIs'].tolist()
# Fetching publication information for each author
pub_details = {}
n = 0
name = ''

for author_name in PI:
    try:
        search_query = scholarly.search_author(author_name)
        author = next(search_query)
        author_info = scholarly.fill(author)

        publications = author_info['publications']
        pub_info = []

        for pub in publications:
            title = pub['bib'].get('title', 'N/A')
            year = pub['bib'].get('pub_year', 'N/A')
            venue = pub['bib'].get('citation', 'N/A')
            name = author_info['name']
            # Gathering publication details for each author
            if re.search(r'\bNature\b', venue) and not re.search(r'\bSpringer\b', venue) and not re.search(r'\bDiscrete\b', venue):
                pub_info.append({'Title': title, 'Author': name, 'Year': year, 'Venue': venue})
        n += 1
    except StopIteration as e:
        # print(f"{PI[n]} not found on Google Scholar.")
        print(e)
        n += 1
        continue

    # Assigning the publication details to the corresponding author
    print(author_name)
    pub_details[author_name] = pub_info

# Writing publication details for each author to a CSV file
file_name = 'MIRMI_full_paperlist.csv'

with open(file_name, mode='w', newline='', encoding='utf-8') as file:
    writer = csv.writer(file)
    writer.writerow(['Publication Title', 'Author', 'Publication Year', 'Venue', 'bib'])

    # Writing publication details separated by author name to the CSV file
    for author, details in pub_details.items():
        for pub in details:
            writer.writerow([pub['Title'], pub['Author'], pub['Year'], pub['Venue']])

print(f"Publication details separated by author names have been written to '{file_name}' successfully.")
