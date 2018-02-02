import pandas
import os
from irp.utils import get_db_summary, get_assemblies, query_reference, filter_hits


def main():

    summary_file = get_db_summary(db="refseq", kingdom="bacteria")

    _ = get_assemblies(species="Staphylococcus aureus", summary=summary_file,
                       outdir="saureus")

    # Path contains both RefSeq and NCTC reference sequences:
    # nctc-tools --species "Staphylococcus aureus" --path ./nctc_saureus make --chromosomes

    files = [os.path.join(os.getcwd(), "saureus", file) for file in os.listdir("saureus") if file.endswith(".fasta")]

    results = []
    for genome in files:
        df, file = query_reference("resources/irp.fasta", reference=genome, outdir="saureus/queries", evalue=1e-06)
        results.append(df)

    df = pandas.concat(results)

    df.to_csv("saureus/hits.tab", sep="\t")

    df.sort_values(by="file", inplace=True)

    df_filtered = filter_hits(df, min_coverage=0.90, min_identity=0.90, min_hits=3)

    df_filtered.sort_values(by="file", inplace=True)

    df_filtered.to_csv("saureus/filtered_hits.tab", sep="\t")


main()