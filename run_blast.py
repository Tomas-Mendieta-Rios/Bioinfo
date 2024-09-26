from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML

def run_blastp(query_file, blast_db, evalue=0.001, outfmt=5, output_file="blast_results.xml"):
    """
    Runs a BLASTP search on the query file against the specified database and parses the results.
    
    Parameters:
    - query_file (str): Path to the FASTA file containing the query sequence(s).
    - blast_db (str): Path to the BLAST database (e.g., Swiss-Prot).
    - evalue (float): E-value threshold for reporting alignments.
    - outfmt (int): Output format (5 is XML).
    - output_file (str): File to save the BLASTP results.
    
    Returns:
    - None: Prints parsed BLAST results.
    """
    
    # Set up the BLASTP command
    blastp_cline = NcbiblastpCommandline(query=query_file, db=blast_db, evalue=evalue, outfmt=outfmt, out=output_file)
    
    # Run the BLASTP command
    print("Running BLASTP search...")
    stdout, stderr = blastp_cline()
    print(f"BLAST search complete. Results saved to {output_file}")
    
    # Parse the BLAST results
    print("Parsing BLAST results...")
    with open(output_file) as result_handle:
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

def blast_fasta(blast_xml_file, output_fasta_file="ordered_blast_output.fasta"):
    """
    Parses BLAST XML results, sorts them by best match (based on e-value),
    and creates a FASTA file from the aligned sequences.
    
    Parameters:
    - blast_xml_file (str): The path to the BLAST results in XML format.
    - output_fasta_file (str): The name of the output FASTA file.
    
    Returns:
    - None: Writes the aligned sequences to a FASTA file ordered by best match.
    """
    results = []
    
    # Step 1: Parse the BLAST results and collect alignments with e-values
    with open(blast_xml_file) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        
        for blast_record in blast_records:
            query_id = blast_record.query
            for alignment in blast_record.alignments:
                for hsp in alignment.hsps:
                    results.append({
                        'title': alignment.title,
                        'evalue': hsp.expect,  # E-value
                        'sequence': hsp.sbjct  # Subject sequence
                    })
    
    # Step 2: Sort results by e-value (from lowest to highest)
    sorted_results = sorted(results, key=lambda x: x['evalue'])
    
    # Step 3: Write sorted results to the output FASTA file
    with open(output_fasta_file, "w") as fasta_file:
        for result in sorted_results:
            fasta_file.write(f">{result['title']} | e-value: {result['evalue']}\n")
            fasta_file.write(f"{result['sequence']}\n")
            fasta_file.write("\n")
    
    print(f"FASTA file '{output_fasta_file}' created successfully, ordered by best match (lowest e-value).")

query_file = "traduccion.fasta"
blast_db = "/usr/local/share/blastdb/swissprot"
run_blastp (query_file,blast_db )
blast_fasta("blast_results.xml", output_fasta_file="ordered_aligned_sequences.fasta")