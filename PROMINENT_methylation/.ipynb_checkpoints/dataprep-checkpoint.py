import argparse
import pandas as pd

def dataprep():
    parser = argparse.ArgumentParser(description='Prepare the data for training. It will generate the sparse matrix for masking.')
    parser.add_argument('--input_csv', default='gene.average.beta.by.intensity.csv', help='Input csv file of X_i. Rows are genes, columns are samples.')
    parser.add_argument('--input_gmt', default='c5.go.bp.v2023.1.Hs.symbols.gmt', help='Input pathway gene sets from MsigDB')
    parser.add_argument('--output', default='pathway_gobp.csv', help='Output file path. csv file.')

    args = parser.parse_args()

    input_file_csv = args.input_csv
    input_file_gmt = args.input_gmt
    output_file = args.output

    # read in X
    m = pd.read_csv(input_file_csv, index_col=0)
    m_genes = list(m.index)
    
    #read in gene sets and get the M
    path_dict = {}
    with open(input_file_gmt, 'r') as file:
        for line in file:
            columns = line.strip().split('\t')
            path = columns[0]
            path_genes = columns[2:]
            if len(list(set(m_genes) & set(path_genes))) >= 10:
                path_dict.update({path:list(set(m_genes) & set(path_genes))})
                
    df = pd.DataFrame(0, index=m_genes, columns = list(path_dict.keys()))
    for path, genes in path_dict.items():
        for gene in genes:
            df.loc[gene,path] = 1
    df = df.reset_index(drop=False)
    df.to_csv(output_file,index=None)

if __name__ == '__main__':
    dataprep()