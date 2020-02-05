import re
import matplotlib.pyplot as plt


class GFF3:
    """GFF3 entry:
    - Hoeveel exonen?
    - Lengte gen
    - Chromosoom
    - Accessiecode
    """

    def __init__(self):
        self.aantal_exonen = 0
        self.lengte_gen = 0
        self.chromosoom = 0
        self.accessiecode = 0

    def set_exonen(self, ex):
        """Bepaal het aantal exonen voor dit gen"""
        self.aantal_exonen = ex

    def get_exonen(self):
        """Return het aantal exonen"""
        return self.aantal_exonen

    def set_lengtegen(self, start, stop):
        """Bepaal de lengte van het gen door stop - start te berekenen

        Input:
        start - str - start van het gen
        stop - str - stop van een gen
        """
        self.lengte_gen = int(stop) - int(start)

    def get_lengtegen(self):
        """Return lengte van het gen"""
        return self.lengte_gen

    def set_chromosoom(self, chr):
        """Bepaal het chromosoom voor dit gen"""
        self.chromosoom = chr

    def get_chromosoom(self):
        """Return het chromosoom van dit gen"""
        return self.chromosoom

    def set_accessiecode(self, accessiecode):
        """Bepaal de accessiecode van dit gen"""
        self.accessiecode = accessiecode

    def get_accessiecode(self):
        """Return de accessiecode van dit gen"""
        return self.accessiecode


def fasta_parser(fasta, regex):
    """Parse FASTA file

    :param fasta: str - name fasta file
    :param regex: str - regular expression kinase
    :return: fasta_dictio - fasta with only sequences with kinase
    expression
    """
    fasta_dictio_temp = {}
    # open and read file
    try:
        with open(fasta) as inFile:
            for line in inFile:
                if line.startswith(">"):
                    key = line.split()[0].split(">")[1]
                    fasta_dictio_temp[key] = []
                else:
                    fasta_dictio_temp[key].append(line.strip())
                    
    except FileNotFoundError:
        print("Dit bestand wordt niet gevonden")
        fasta = input("Geef een nieuw bestand: ")
        fasta_parser(fasta, regex)

    # get only kinase sequences instead of all of them
    fasta_dictio = get_kinase_sequences(fasta_dictio_temp, regex)
    return fasta_dictio


def get_kinase_sequences(fasta_dictio_temp, regex):
    """Select de sequenties met deze kinaze regex

    :param fasta_dictio_temp: dict - temporarily dictionary with
    sequences
    :param regex: str - regular expression of kinase
    :return: fasta_dictio - dict - dictionary with only sequences with
    kinase
    """
    fasta_dictio = {}
    for key, value in fasta_dictio_temp.items():
        if find_kinase("".join(value), regex):
            fasta_dictio[key] = "".join(value)
    return fasta_dictio


def find_kinase(sequence, regex):
    """Find kinase in sequence

    :param sequence: sequence to look for kinase
    :param regex: regular expression
    :return: True/False
    """
    match = re.search(regex, sequence)
    if match:
        return True
    else:
        return False


def gff3_parser(gff3, fasta_dic):
    """Create a list with objects that match all fasta entries plus
    additional info like:
        - Hoeveel exonen?
        - Lengte gen
        - Chromosoom
        - Accessiecode
    :param gff3: str - naam bestand
    :param fasta_dic: dict - fasta dictionary {header: seq}
    :return: list with gff3 objects
    """
    gff3_entries = []
    entry = ""
    # open en lees bestand
    try: 
        with open(gff3) as inFile:
            for line in inFile:
                # als mRNA in de regel staat
                if re.search("mRNA", line):
                    # als entry niet leeg is
                    if entry != "":
                        entry.set_exonen(exons)
                        if entry.get_accessiecode() in fasta_dic.keys():
                            gff3_entries.append(entry)
                            exons = 0
                        # Maak een nieuwe entry aan
                        entry = GFF3()
                        entry.set_chromosoom(line.split()[0].split("Chr")[1])
                        entry.set_lengtegen(line.split()[3], line.split()[4])
                        entry.set_accessiecode(line.split()[8].split(";")[0] \
                                           .split("=")[1])
                    # als entry leeg is, maak een nieuwe aan
                    else:
                        exons = 0
                        entry = GFF3()
                        # Haal chromosoom, start, stop en accessiecode uit
                        # de regel en sla deze op in het object
                        entry.set_chromosoom(line.split()[0].split("Chr")[1])
                        entry.set_lengtegen(line.split()[3], line.split()[4])
                        entry.set_accessiecode(line.split()[8].split(";")[0]\
                                           .split("=")[1])
                # Tel het aantal exonen
                elif re.search("exon", line):
                    exons += 1
        except FileNotFoundError:
            print("Dit bestand wordt niet gevonden")
            gff3 = input("Geef een nieuw bestand: ")
            gff3_parser(gff3, fasta_dic)

    return gff3_entries


def matplotlib_data_prep(gff3_ol):
    """Prep matplotlib graph data zodat deze in een grafiek weergegeven
    kan worden

    :param gff3_ol: list with gff3 objects
    :return: matplotlib_data: dict - {chr1: count, chr2: count...}
    """
    matplotlib_data = {}
    # voor ieder object in de lijst
    for obj in gff3_ol:
        # sla het chromosoomnummer op
        chr = obj.get_chromosoom()
        # als dit chromosoom niet in matplotlib_data zit, voegen we
        # deze toe, plus een counter (staat op 1)
        if chr not in matplotlib_data:
            matplotlib_data[chr] = 1
        # als dit chromosoom wel in matplotlib_data zit, tellen we een
        # op bij de counter
        else:
            matplotlib_data[chr] += 1

    return matplotlib_data


def create_graph(matplotlib_data):
    """Create matplotlib graph

    :param matplotlib_data:
    :return: /
    """
    # Maak de barplot
    plt.bar(list(matplotlib_data.keys()), list(matplotlib_data.values()))
    plt.xlabel("Chromosome number")
    plt.ylabel("Counts")
    plt.title("Gene counts for each chromosome with kinase domain")
    plt.show()
    return


def output_no_gui(gff3_ol):
    """Pretty print van de output zonder GUI

    :param gff3_ol:
    :return:
    """
    # zet doorgaan op True
    doorgaan = True
    # print iedere accessiecode
    for el in gff3_ol:
        print(el.get_accessiecode())
    # als doorgaan gelijk is aan True blijft deze loop doorgaan
    while doorgaan:
        accession = input("Kies een accessiecode")
        for el in gff3_ol:
            if accession == el.get_accessiecode():
                print("Aantal exonen: ", el.get_exonen())
                print("Lengte gen: ", el.get_lengtegen())
                print("Chromosoom: ", el.get_chromosoom())
        vraag = input("Wilt u nog doorgaan? ")
        if vraag == "Ja" or vraag == "ja":
            doorgaan = True
        else:
            doorgaan = False


if __name__ == "__main__":
    # declareer variabelen
    gff3 = "TAIR10_GFF3_genes.gff"
    fasta = "TAIR10_pep_20101214.fa"
    regex = "[LIVMFYC].[HY].D[LIVMFY]K.{2}N[LIVMFYCT]{3}"

    # Lees de bestanden in
    fasta_filtered_dic = fasta_parser(fasta, regex)
    gff3_ol = gff3_parser(gff3, fasta_filtered_dic)

    # Creeër de grafiek + dataprep
    matplotlib_data = matplotlib_data_prep(gff3_ol)
    create_graph(matplotlib_data)

    # Print de output zonder GUI
    output_no_gui(gff3_ol)
