from scholarly import scholarly
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

nature_keywords = [
    'Nature communications',
    'Nature Machine Intelligence',
    'Nature Biomedical Engineering',
    'Nature Neuroscience',
    'Nature Biomedical Engineering',
]

keywords_alter = [
    'Scientific Reports',#Nature
    'Scientific data',#Nature
    'Communications Engineering',#Nature
    'Microsystems & Nanoengineering',#Nature
    'npj Robotics',#Nature
    'npj Science of Learning',#Nature
    'Cyborg and Bionic Systems'#Science
]

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
            #r'\bNature\s*\d+\b'
            if re.search(r'\bNature\b', venue) and not re.search(r'\bSpringer\b', venue) and not re.search(r'\bDiscrete\b', venue): #To search only Nature related papers 
                if any(re.search(fr'\b{keyword}\b', venue, re.IGNORECASE) for keyword in nature_keywords):
                    pub_info.append({'Author': name, 'Year': year, 'Venue': venue, 'Title': title, 'Affiliation' : 'TUM'})
            if any(re.search(fr'\b{keyword}\b', venue, re.IGNORECASE) for keyword in keywords_alter):
                pub_info.append({'Author': name, 'Year': year, 'Venue': venue, 'Title': title, 'Affiliation' : 'TUM'})
            if re.search(r'\bNature\s*\d+\b', venue):
                pub_info.append({'Author': name, 'Year': year, 'Venue': venue, 'Title': title, 'Affiliation' : 'TUM'})
            # if re.search(r'\bNature\b', venue) and not re.search(r'\bSpringer\b', venue) and not re.search(r'\bDiscrete\b', venue): #To search only Nature related papers
        n += 1
        pub_details[author_name] = pub_info
    except StopIteration as e:
        # print(f"{PI[n]} not found on Google Scholar.")
        print(e)
        n += 1
        continue

    # Assigning the publication details to the corresponding author
    print(author_name)

# Creating a DataFrame from the gathered publication details
all_pub_details = []

for author_name, pub_info in pub_details.items():
    all_pub_details.extend(pub_info)

df_pub_details = pd.DataFrame(all_pub_details)

# Save the DataFrame to an Excel file
output_file_path = './publication_details.xlsx'  # Replace with your desired output file path
df_pub_details.to_excel(output_file_path, index=False)

print(f"Publication details saved to {output_file_path}")
# Writing publication details for each author to an xlsx file
# file_name = 'MIRMI_full_paperlist.xlsx'

# with pd.ExcelWriter(file_name, engine='xlsxwriter') as writer:
#     for author, details in pub_details.items():
#         df = pd.DataFrame(details)
#         # Adding 'TUM' column with constant value 'TUM'
#         df['Affiliation'] = 'TUM'
#         df.to_excel(writer, sheet_name=author, index=False)

# print(f"Publication details separated by author names have been written to '{file_name}' successfully.")
# # Writing publication details for each author to a CSV file
# file_name = 'MIRMI_full_paperlist.csv'

# with open(file_name, mode='w', newline='', encoding='utf-8') as file:
#     writer = csv.writer(file)
#     writer.writerow(['Publication Title', 'Author', 'Publication Year', 'Venue', 'Affiliation'])

#     # Writing publication details separated by author name to the CSV file
#     for author, details in pub_details.items():
#         for pub in details:
#             writer.writerow([pub['Title'], pub['Author'], pub['Year'], pub['Venue'], 'TUM'])

# print(f"Publication details separated by author names have been written to '{file_name}' successfully.")
