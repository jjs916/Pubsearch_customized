import dblp
import time
import pandas as pd

max_retries = 3  # 최대 시도 횟수
search_timeout = 3  # 검색 타임아웃 (초)
buffer_time = 1.2

# Read Excel file and extract author names
file_path = './PIlist.xlsx'  # Replace with your file path
df = pd.read_excel(file_path)

# Assuming the column containing author names is named 'Authors'
PI = df['PIs'].tolist()
print(f"Total number of PIs: {str(len(PI))}")
print(f"Estimate Time: {str((len(PI)*(1))/3600)} h")
# publication saving list
publications_data = []

programstarttime = time.time()
totpub = 0
for authorname in PI: #repeat for all PIs
    pubdata = []
    n=0
    retries = 0
    while retries < max_retries:
    #do a simple author search
        try:
            start_time = time.time()
            authors = dblp.search(authorname)
            elapsed_time = time.time() - start_time
            if elapsed_time > search_timeout:
                print(f"Search timeout reached ({search_timeout} seconds). Stopping the process.")
                break
            if authors:
                author = authors[0]
                author.load_data
                print(author.name)
                # author.publications[0].load_data()
                # max = range(len(author.publications))
                # print(authorname, len(author.publications))
                # totpub += len(author.publications)
                publications_data.append({'Author': author.name})
                break
            else:
                break
        except Exception as e:
            print(f"Error occurred for {authorname}: {str(e)}")
            retries += 1
            time.sleep(60)
            
    time.sleep(1.5)

# make dataframe
df = pd.DataFrame(publications_data)

# save excel file
df.to_excel(f'PIname.xlsx', index=False)
programendtime = time.time() - programstarttime
print('Searching Finished!')
print('Total '+str(programendtime)+' spent')
print('Total '+str(totpub)+' publications')