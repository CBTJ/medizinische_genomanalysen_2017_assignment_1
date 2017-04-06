import mysql.connector
import pysam
import pybedtools

# fuer die Durchfuehrung braucht man Internetzugang um auf die DB zugreifen zu koennen
# python 3
# das BAM-File
# das indexierte BAM-file = .bai (samtools)

__author__ = 'Claudia Juno'

class Genedetails:                  # Klasse um Daten der Internetabfrage (fetch_gene_coordinates) zu speichern
    def __init__(self):
        self.genename2 = ""
        self.genename = ""
        self.genechrom = ""
        self.genetxstart = 0
        self.genetxend = 0
        self.genestrand = ""
        self.exoncount = 0
        self.exonstart = ""
        self.exonends = ""

class Assignment1:
    
    def __init__(self,f):
        ## Your gene of interest

        self.gene = "WT1"
        self.genedetails = Genedetails()
        self.fetch_gene_coordinates("hg19", "assignment1fetch.TXT")


        ## bam file

        self.bamfile = f
        self.samfile = pysam.AlignmentFile(f, "rb") # aus BAM wird SAM
        self.lines = list(self.samfile.fetch(self.genedetails.genechrom, self.genedetails.genetxstart,
                                             self.genedetails.genetxend)) # Abfrage

        #fuer Berechnung der beiden coverages
        a = pybedtools.BedTool(self.bamfile)
        self.b = a.genome_coverage(bg=True)


    def fetch_gene_coordinates(self, genome_reference, file_name):
        
        print ("Connecting to UCSC to fetch data")
        
        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=genome_reference)
        
        ## Get cursor
        cursor = cnx.cursor()
        
        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        ## Build query

        query = "SELECT DISTINCT %s from refGene" % ",".join(query_fields) + \
                " WHERE refGene.name2=" + '"' + self.gene + '"' ""
        
        ## Execute query
        cursor.execute(query)
        
        ## Write to file
        ## Speichern der Daten aus der DB Abfrage in der Klasse genedetails
        with open(file_name, "w") as fh:
            for row in cursor:
                if row [0] == self.gene:
                    fh.write(str(row) + "\n")

                    self.genedetails.genename2 = row[0]
                    self.genedetails.genename = row[1]
                    self.genedetails.genechrom = row[2][3:]
                    self.genedetails.genetxstart = row[3]
                    self.genedetails.genetxend = row[4]
                    self.genedetails.genestrand = row[5]
                    self.genedetails.exoncount = row[6]
                    self.genedetails.exonstart = row[7]
                    self.genedetails.exonends = row[8]
            
        ## Close cursor & connection
        cursor.close()
        cnx.close()
        
        print ("Done fetching data")
                
    def get_sam_header(self):

        for a, b  in self.samfile.header['HD'].items():
            print(a, b)
            #optional fuer SQ als Teil des headers, einkommentieren
        #for a in self.samfile.header ['SQ']:
            #print (a)


    def get_properly_paired_reads_of_gene(self):
        c = 0
        z = []
        for read in self.lines:
            if read.is_paired:
                z.append (read)
                c += 1
        print (c)
        # wenn nicht Anzahl sondern Zeilen erwuenscht bitte die beiden folgenden Zeilen einkommentieren
        #for i in z:
            #print (i)


    def get_gene_reads_with_indels(self):
        c = 0
        cigarlist = []
        for l in self.lines:
            if (not (l.is_unmapped)): #ergibt die mapped reads
                cigarLine = l.cigar
                for (cigarType, cigarLine) in cigarLine:
                    if (cigarType == 1) or (cigarType == 2): # cigar Type 1 = insertion, 2 = deletion
                        cigarlist.append (l)
                        c += 1
        print (c)
        # wenn nicht Anzahl sondern Zeilen erwuenscht bitte die beiden folgenden Zeilen einkommentieren
        #for c in cigarlist:
            #print (c)


    def calculate_total_average_coverage(self):

        z = 0
        sum = 0

        for line in self.b:                 # self.b genome coverage Berechnung ueber BedTools
            sum += float(line[3])
            z += 1
            #print (line)

        tacoverage = sum / z
        print (round (tacoverage, 4))

    def calculate_gene_average_coverage(self):

        z = 0
        sum = 0
        gacoverage = 0

        for line in self.b:     # self.b genome coverage Berechnung ueber BedTools
            c = int (line[1])
            d = int (line [2])
            if self.genedetails.genetxstart <= c and self.genedetails.genetxend >= d: # eingrenzen da gene coverage
                sum += float(line[3])
                z += 1
        if z > 0:
            gacoverage = sum / z
        print(round(gacoverage, 4))

        
    def get_number_mapped_reads(self):
        z = 0
        for l in self.lines:
            if (not (l.is_unmapped)): #ergibt die mapped reads
                z +=1
        print(z)
        
    def get_gene_symbol(self):
        print(self.genedetails.genename2)
        
    def get_region_of_gene(self):
        print("Chromosome {} from {} to {}".format (self.genedetails.genechrom, self.genedetails.genetxstart,
                                                    self.genedetails.genetxend))
        
    def get_number_of_exons(self):
        print(self.genedetails.exoncount)


    def print_summary(self):
        print ("-----------------------")
        print ("Header:")
        self.get_sam_header()
        print("-----------------------")
        print("Properly paired reads:")
        self.get_properly_paired_reads_of_gene()
        print("-----------------------")
        print("Reads with indels:")
        self.get_gene_reads_with_indels()
        print("-----------------------")
        print("Total average coverage:")
        self.calculate_total_average_coverage()
        print("-----------------------")
        print("Gene average coverage:")
        self.calculate_gene_average_coverage()
        print("-----------------------")
        print("Number of mapped reads:")
        self.get_number_mapped_reads()
        print("-----------------------")
        print("Gene symbol:")
        self.get_gene_symbol()
        print("-----------------------")
        print("Region of Gene:")
        self.get_region_of_gene()
        print("-----------------------")
        print("Number of exons:")
        self.get_number_of_exons()
        print("-----------------------")

if __name__ == '__main__':
    print ("Assignment 1")
    print(__author__)
    assignment1 = Assignment1('HG00096.chrom11.ILLUMINA.bwa.GBR.low_coverage.20120522.bam')
    assignment1.print_summary()
    

    
    

