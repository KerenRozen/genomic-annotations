import pandas as pd
import  time


class GeneDetector:

    # input is the path to the genome DB
    def __init__(self, path):
        self.genes_db = self.process_data(path)
        #self.n = self.genes_db.shape[0]
        self.chromosomes_range = self.map_chromosomes_rows({})

    def process_data(self, path):
        cols = ["seqid", "type", "start", "end"]
        df = pd.read_csv(path, usecols=[0, 2, 3, 4],  names=cols, dtype={'seqid': 'str'}, comment='#', sep='\t', header=None)
        df = df[df['type'] == 'gene'].reset_index(drop=True)
        #df = df.replace(to_replace=r'NC_0{4,5}(\d{1,2})\.\d+', value=r'\1', regex=True)

        return df

    def map_chromosomes_rows(self, d):
        chromosome_to_int = {'10': 10, '11': 11, '12': 12, '13': 13, '14': 14, '15': 15, '16': 16, '17': 17, '18': 18,
                             '19': 19, '1': 1, '20': 20, '21': 21, '22': 22, '2': 2, '3': 3, '4': 4, '5': 5, '6': 6,
                             '7': 7, '8': 8, '9': 9, 'MT': 23, 'X': 24, 'Y': 25, "GL000192.1": 26, "GL000225.1": 27,
                             "GL000194.1": 28, "GL000193.1": 29, "GL000200.1": 30, "GL000222.1": 31, "GL000212.1": 32,
                             "GL000195.1": 33, "GL000223.1": 34, "GL000224.1": 35, "GL000219.1": 36, "GL000205.1": 37,
                             "GL000215.1": 38, "GL000216.1": 39, "GL000217.1": 40, "GL000199.1": 41, "GL000211.1": 42,
                             "GL000213.1": 43, "GL000220.1": 44, "GL000218.1": 45, "GL000209.1": 46, "GL000221.1": 47,
                             "GL000214.1": 48, "GL000228.1": 49, "GL000227.1": 50, "GL000191.1": 51, "GL000208.1": 52,
                             "GL000198.1": 53, "GL000204.1": 54, "GL000233.1": 55, "GL000237.1": 56, "GL000230.1": 57,
                             "GL000242.1": 58, "GL000243.1": 59, "GL000241.1": 60, "GL000236.1": 61, "GL000240.1": 62,
                             "GL000206.1": 63, "GL000232.1": 64, "GL000234.1": 65, "GL000202.1": 66, "GL000238.1": 67,
                             "GL000244.1": 68, "GL000248.1": 69, "GL000196.1": 70, "GL000249.1": 71, "GL000246.1": 72,
                             "GL000203.1": 73, "GL000197.1": 74, "GL000245.1": 75, "GL000247.1": 76, "GL000201.1": 77,
                             "GL000235.1": 78, "GL000239.1": 79, "GL000210.1": 80, "GL000231.1": 81, "GL000229.1": 82,
                             "GL000226.1": 83, "GL000207.1": 84}

        for char in chromosome_to_int.keys():
            rows = self.genes_db[self.genes_db['seqid'].astype(str) == char].index.tolist()
            if rows:
                d[chromosome_to_int[char]] = (rows[0], rows[-1])
            else:
                d[chromosome_to_int[char]] = (0, 0)

        return d

    def inside_gene(self, chromosome, start_pos, end_pos):
        if end_pos < start_pos or start_pos < 1:
            return False

        left = self.chromosomes_range[chromosome][0]
        right = self.chromosomes_range[chromosome][1]

        while left <= right:
            mid = (left+right) // 2
            if self.genes_db.iloc[mid, 2] <= start_pos and end_pos <= self.genes_db.iloc[mid, 3]:
                return True
            if self.genes_db.iloc[mid, 2] > start_pos:
                right = mid - 1
            else:
                left = mid + 1

        return False



db = GeneDetector('C:/Users/keren/PycharmProjects/project/Homo_sapiens.GRCh37.87.gff3')
data = pd.read_csv('../data/example_read_file.tsv', usecols=[1, 2, 3], sep='\t')
#data['length'] = data['end_pos'] - data['start_pos']
#print(data.length.max())
#data.to_csv("data", sep='\t', encoding='utf-8')


def tester(df):
    for index, row in df.iterrows():
        if row['chromosome'] <= 84:
            print(index)
            print(db.inside_gene(row['chromosome'], row['start_pos'], row['end_pos']))

test_dict = {i:i for i in range(200)}



start_time = time.time()

for i in range(200):
    print(test_dict[i])


print("--- %s seconds ---" % (time.time() - start_time))