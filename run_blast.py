from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

# Step 1: Define query file and database paths
query_file = "traduccion.fasta"  # Replace with your sequence file
blast_db = "/usr/local/share/blastdb/swissprot"  # Path to the Swiss-Prot database

# Step 2: Set up the BLASTP command
blastp_cline = NcbiblastpCommandline(query=query_file, db=blast_db, evalue=0.001, outfmt=5, out="blast_results.xml")

# Step 3: Run the BLASTP command
print("Running BLASTP search...")
stdout, stderr = blastp_cline()
print("BLAST search complete. Results saved to blast_results.xml")

# Step 4: Parse the BLAST results
print("Parsing BLAST results...")
with open("blast_results.xml") as result_handle:
    blast_records = NCBIXML.parse(result_handle)
    
    for blast_record in blast_records:
        print(f"Query: {blast_record.query}")
        for alignment in blast_record.alignments:
            print("****Alignment****")
            print(f"Sequence: {alignment.title}")
            print(f"Length: {alignment.length}")
            for hsp in alignment.hsps:
                print(f"e-value: {hsp.expect}")
                print(f"Query: {hsp.query}")
                print(f"Match: {hsp.match}")
                print(f"Subject: {hsp.sbjct}")
                print()

print("Parsing complete.")
