import json
import requests
import numpy as np
import pandas as pd
# temporay gene list

def Get_Gene_array(filename) :
    
    '''filename = "example.txt"'''
    
    genes = np.loadtxt(filename, dtype='object', unpack=True)
    return genes
# temporay gene list

def Get_background_array(filename) :
    
    '''filename = "example.txt"'''
    
    background = np.loadtxt(filename, dtype='object', unpack=True)
    return background
# Analyze gene set

def write_enrichment_file(genes, database, export_filename) :
    
    '''example: database = "KEGG_2021_Human"
       example: export_filename = "test"'''
    
    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/addList'
    genes_str = '\n'.join(genes)
    description = 'Example gene list'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    response = requests.post(ENRICHR_URL, files=payload)
    if not response.ok:
        raise Exception('Error analyzing gene list')

    data = json.loads(response.text)
    print(data)

    ENRICHR_URL = 'https://maayanlab.cloud/Enrichr/export'
    query_string = '?userListId=%s&filename=%s&backgroundType=%s'
    user_list_id = data['userListId']
    filename = export_filename
    gene_set_library = database

    url = ENRICHR_URL + query_string % (user_list_id, filename, gene_set_library)
    response = requests.get(url, stream=True)
    print(url)
    with open(filename + '.txt', 'wb') as f:
        for chunk in response.iter_content(chunk_size=1024): 
            if chunk:
                f.write(chunk)

# Analyze gene set

def write__backgroud_enrichment_file(genes,background, database, export_filename) :
    

    base_url = "https://maayanlab.cloud/speedrichr"

    description = "sample gene set with background"

    res = requests.post(
        base_url+'/api/addList',
        files=dict(
          list=(None, '\n'.join(genes)),
          description=(None, description),
        )
      )
    if res.ok:
        userlist_response = res.json()
        print(userlist_response)
        
    
    base_url = "https://maayanlab.cloud/speedrichr"



    res = requests.post(
        base_url+'/api/addbackground',
        data=dict(background='\n'.join(background)),
      )

    if res.ok:
        background_response = res.json()
        print(background_response)
    
    base_url = "https://maayanlab.cloud/speedrichr"
    gene_set_library = database

    res = requests.post(
            base_url+'/api/backgroundenrich',
            data=dict(
            userListId=userlist_response['userListId'],
            backgroundid=background_response['backgroundid'],
            backgroundType=gene_set_library,
            )
        )
    if res.ok:
        results = res.json()
    #    print(results)
    df = pd.DataFrame(results)

    df['Term'] = df[gene_set_library].apply(lambda x: x[1])
    df['P-value'] = df[gene_set_library].apply(lambda x: x[2])
    df['Genes'] = df[gene_set_library].apply(lambda x: ';'.join(x[5]))
    df['Adjusted P-value'] = df[gene_set_library].apply(lambda x: x[6])

    df.drop(columns=[gene_set_library], inplace=True)

    # Define the file path where you want to save the CSV
    filename = export_filename
    
    df.to_csv(filename + '.txt', sep='\t', index=False)
                                            

import requests
import json
import pandas as pd
import io

def get_enrichment_dataframe(genes, database):
    """
    Analyze a gene set using Enrichr and return the enrichment results as a pandas DataFrame.

    Parameters:
    - genes (list of str): List of gene symbols.
    - database (str): Name of the Enrichr database to use (e.g., "KEGG_2021_Human").

    Returns:
    - pd.DataFrame: DataFrame containing the enrichment results.
    
    Example:
        enrichment_df = get_enrichment_dataframe(
            genes=["BRCA1", "TP53", "EGFR"],
            database="KEGG_2021_Human"
        )
    """
    
    # Endpoint URLs
    ADDLIST_URL = 'https://maayanlab.cloud/Enrichr/addList'
    EXPORT_URL = 'https://maayanlab.cloud/Enrichr/export'
    
    # Prepare the gene list payload
    genes_str = '\n'.join(genes)
    description = 'Gene list for enrichment analysis'
    payload = {
        'list': (None, genes_str),
        'description': (None, description)
    }

    # Submit the gene list to Enrichr
    response = requests.post(ADDLIST_URL, files=payload)
    if not response.ok:
        raise Exception(f'Error adding gene list to Enrichr: {response.text}')

    data = response.json()
    print("Enrichr Add List Response:", data)

    user_list_id = data.get('userListId')
    if not user_list_id:
        raise Exception('No userListId found in Enrichr response.')

    # Prepare the export request
    query_string = f'?userListId={user_list_id}&filename=enrichment_results&backgroundType={database}'
    export_full_url = EXPORT_URL + query_string
    print("Export URL:", export_full_url)

    # Request the enrichment results
    response = requests.get(export_full_url, stream=True)
    if not response.ok:
        raise Exception(f'Error fetching enrichment results: {response.text}')

    # Read the response content
    content = response.content.decode('utf-8')

    # Use StringIO to read the content into pandas as if it were a file
    df = pd.read_csv(io.StringIO(content), sep='\t')

    return df
                


def write_background_enrichment_df(genes, background, database):
    base_url = "https://maayanlab.cloud/speedrichr"
    description = "sample gene set with background"

    # Add gene list
    res = requests.post(
        f"{base_url}/api/addList",
        files={
            'list': (None, '\n'.join(genes)),
            'description': (None, description),
        }
    )
    if not res.ok:
        raise Exception(f"Failed to add gene list: {res.text}")
    
    userlist_response = res.json()
    print("User List Response:", userlist_response)

    # Add background
    res = requests.post(
        f"{base_url}/api/addbackground",
        data={'background': '\n'.join(background)}
    )
    if not res.ok:
        raise Exception(f"Failed to add background: {res.text}")
    
    background_response = res.json()
    print("Background Response:", background_response)

    # Perform background enrichment
    res = requests.post(
        f"{base_url}/api/backgroundenrich",
        data={
            'userListId': userlist_response['userListId'],
            'backgroundid': background_response['backgroundid'],
            'backgroundType': database,
        }
    )
    if not res.ok:
        raise Exception(f"Failed to perform background enrichment: {res.text}")
    
    results = res.json()
    
    # Create DataFrame
    df = pd.DataFrame(results)
    df['Term'] = df[database].apply(lambda x: x[1])
    df['P-value'] = df[database].apply(lambda x: x[2])
    df['Genes'] = df[database].apply(lambda x: ';'.join(x[5]))
    df['Adjusted P-value'] = df[database].apply(lambda x: x[6])

    # Remove the original column
    df.drop(columns=[database], inplace=True)

    return df


