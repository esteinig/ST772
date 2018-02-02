from Bio import Entrez
from urllib.request import urlretrieve
from urllib.error import HTTPError, URLError

import os
import pandas
import gzip
import shutil

from pandas.errors import EmptyDataError

import subprocess

Entrez.email = "pathfinder@jcu.edu.au"


def get_taxid(species):

    """https://gist.github.com/brantfaircloth/802094"""

    species = species.replace(" ", "+").strip()
    search = Entrez.esearch(term=species, db="taxonomy", retmode="xml")
    record = Entrez.read(search)

    try:
        return int(record['IdList'][0])
    except KeyError or ValueError:
        print("Could not find valid TaxID for", species)
        return None


def get_db_summary(db="refseq", kingdom="bacteria", force=False):

    base = os.path.abspath(__file__ + "/../../")
    resources = os.path.join(base, "resources")

    summary = "ftp://ftp.ncbi.nlm.nih.gov/genomes/{}/{}/assembly_summary.txt".format(db, kingdom)

    summary_file = os.path.join(resources, "{}_{}.txt".format(db, kingdom))

    if os.path.exists(summary_file) and not force:
        print("Summary file for {} in {} exists at {}".format(kingdom, db, summary_file))
    else:
        urlretrieve(summary, summary_file)

    return summary_file


def get_assemblies(species, summary, outdir=os.getcwd(), tax_id=None, complete=True):

    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    if not tax_id:
        tax_id = get_taxid(species=species)

    df = pandas.read_csv(summary, sep="\t", header=1)
    df = df.rename(columns={"# assembly_accession": "assembly_accession"})

    if complete:
        df = df[(df["species_taxid"] == tax_id) & (df["assembly_level"] == "Complete Genome")]
    else:
        df = df[(df["species_taxid"] == tax_id)]

    print("Found {} genome assemblies for species {}".format(len(df), species))

    files = []
    for i, row in df.iterrows():
        ftp_path = row["ftp_path"]
        ftp_name = os.path.basename(row["ftp_path"])
        ftp_link = ftp_path + "/" + ftp_name + "_genomic.fna.gz"

        accession = row["assembly_accession"]

        assembly_file = os.path.join(outdir, accession)

        if not os.path.exists(assembly_file+".fasta"):
            print("Downloading", assembly_file+".fna.gz")
            try:
                urlretrieve(ftp_link, assembly_file+".fna.gz")
                print("Decompressing to", assembly_file + ".fasta")
                with gzip.open(assembly_file + ".fna.gz", 'rb') as file_in, \
                        open(assembly_file + ".fasta", 'wb') as file_out:
                    shutil.copyfileobj(file_in, file_out)

                os.remove(assembly_file + ".fna.gz")
            except URLError or HTTPError:
                print("Failed to download", assembly_file)
                continue

        files.append(assembly_file+".fasta")

    return files


def query_reference(query, reference, outdir=os.getcwd(), evalue=1e-06, force=False):

    outdir = os.path.abspath(outdir)
    os.makedirs(outdir, exist_ok=True)

    name, ext = os.path.splitext(os.path.basename(reference))

    db = os.path.join(outdir, name+".db")
    results = os.path.join(outdir, name+".results")

    db_command = 'makeblastdb -in {} -dbtype nucl -out {}'.format(reference, db)

    query_command = 'blastn -query {} -evalue {} -db {} -outfmt "6 qacc sacc pident qlen length ' \
                    'evalue qstart qend sstart send" -out {}'.format(query, evalue, db, results)

    if os.path.exists(results) and not force:
        print("Results from query exist at", results)
        pass
    else:
        try:
            subprocess.call(db_command, shell=True)
        except subprocess.CalledProcessError:
            print("Failed to run makeblastdb with command", db_command)
        except OSError:
            print("Could not detect makeblastdb in $PATH")

        try:
            subprocess.call(query_command, shell=True)

            for ext in (".nhr", ".nin", ".nsq"):
                os.remove(db + ext)

        except subprocess.CalledProcessError:
            print("Failed to run blastn with command", db_command)

    try:
        df = pandas.read_csv(results, sep="\t", header=None)

        df.columns = ["query", "reference", "identity", "query_length", "ref_length", "evalue",
                      "qstart", "qend", "rstart", "rend"]

        df["coverage"] = df["ref_length"] / df["query_length"]
        df["file"] = [os.path.basename(reference) for i in range(len(df))]
    except EmptyDataError:
        df = pandas.DataFrame()

    return df, results


def filter_hits(df, min_coverage=0.90, min_identity=0.90, min_hits=2):

    df = df.loc[(df["coverage"] >= min_coverage) & (df["identity"] >= min_identity*100)]

    if min_hits:
        counts = df["file"].value_counts()
        df = df.loc[df["file"].isin(counts[counts >= min_hits].index.values)]

    return df


