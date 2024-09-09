import ssl
import signal
import subprocess
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Disable SSL certificate verification (for development purposes only; not recommended in production)
ssl._create_default_https_context = ssl._create_unverified_context

# Define a timeout handler to raise an exception if the BLAST search takes too long
def timeout_handler(signum, frame):
    raise TimeoutError("The BLAST search took too long and was terminated.")

# Set the timeout duration (e.g., 300 seconds)
signal.signal(signal.SIGALRM, timeout_handler)

# Function for remote BLAST search
def run_blast_remote(fasta_sequence):
    print("Starting remote BLAST search...")
    signal.alarm(300)  # Set an alarm for 300 seconds
    result_handle = NCBIWWW.qblast("blastp", "nr", fasta_sequence)
    signal.alarm(0)  # Disable the alarm after successful completion
    print("Remote BLAST search completed.")

    # Save the BLAST results to an XML file
    output_file = "blast_remote_result.xml"
    print(f"Saving BLAST results to '{output_file}'...")
    with open(output_file, "w") as out_handle:
        out_handle.write(result_handle.read())
    print("BLAST results saved successfully.")
    result_handle.close()

    return output_file

# Function for local BLAST search
def run_blast_local(fasta_file, db_name="my_local_db"):
    output_file = "blast_local_result.xml"
    blast_command = [
        "blastp",  # Program name (blastp for protein, blastn for nucleotide)
        "-query", fasta_file,  # Input sequence file
        "-db", db_name,  # Local BLAST database name
        "-out", output_file,  # Output file
        "-outfmt", "5",  # Output format (5 for XML)
    ]

    print("Running BLAST locally...")
    subprocess.run(blast_command, check=True)
    print(f"Local BLAST search completed. Results saved to {output_file}.")

    return output_file

# Function to parse BLAST results
def parse_blast_results(xml_file):
    print("Parsing BLAST result XML file...")
    with open(xml_file) as result_handle:
        blast_record = NCBIXML.read(result_handle)

    # Print basic information about the first alignment
    print("Printing alignments...")
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            print("****Alignment****")
            print(f"sequence: {alignment.title}")
            print(f"length: {alignment.length}")
            print(f"e-value: {hsp.expect}")
            print(f"score: {hsp.score}")

# Main script
def main():
    # Sequence in FASTA format (replace it with your sequence if needed)
    fasta_sequence = """>sp|P01308|INS_HUMAN Insulin precursor OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1
    MALWMRLLPLLALLALWGPDPAAAFFVQPCSQITTCGILICSSLSTNQVEALYLVCGERGFFYTPKT
    RRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"""

    # Save the FASTA sequence to a file for local BLAST
    fasta_file = "input_sequence.fasta"
    with open(fasta_file, "w") as f:
        f.write(fasta_sequence)

    # Choose between local and remote BLAST
    print("Choose BLAST method:")
    print("1: Remote BLAST (NCBI)")
    print("2: Local BLAST")
    choice = input("Enter your choice (1 or 2): ")

    if choice == "1":
        # Run remote BLAST search
        xml_file = run_blast_remote(fasta_sequence)
    elif choice == "2":
        # Run local BLAST search
        xml_file = run_blast_local(fasta_file)
    else:
        print("Invalid choice. Please enter 1 or 2.")
        return

    # Parse and display the results
    parse_blast_results(xml_file)

if __name__ == "__main__":
    main()
