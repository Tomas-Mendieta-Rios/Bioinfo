import ssl
import signal
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# Disable SSL certificate verification (for development purposes only; not recommended in production)
ssl._create_default_https_context = ssl._create_unverified_context

# Sequence in FASTA format (feel free to replace it with your original sequence)
fasta_sequence = """>sp|P01308|INS_HUMAN Insulin precursor OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1
MALWMRLLPLLALLALWGPDPAAAFFVQPCSQITTCGILICSSLSTNQVEALYLVCGERGFFYTPKT
RRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"""

# Define a timeout handler to raise an exception if the BLAST search takes too long
def timeout_handler(signum, frame):
    raise TimeoutError("The BLAST search took too long and was terminated.")

# Set the timeout duration (e.g., 300 seconds)
signal.signal(signal.SIGALRM, timeout_handler)
signal.alarm(300)  # Set an alarm for 300 seconds

# Start the BLAST search
print("Starting BLAST search...")
result_handle = NCBIWWW.qblast("blastp", "nr", fasta_sequence)
signal.alarm(0)  # Disable the alarm after successful completion
print("BLAST search completed.")

# Save the BLAST results to an XML file
print("Saving BLAST results to 'blast_remote_result.xml'...")
with open("blast_remote_result.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
print("BLAST results saved successfully.")

# Close the result handle
result_handle.close()
print("Result handle closed.")

# Read and parse the XML file with the BLAST results
print("Parsing BLAST result XML file...")
with open("blast_remote_result.xml") as result_handle:
    blast_record = NCBIXML.read(result_handle)
print("Parsing completed.")

# Print basic information about the first alignment
print("Printing alignments...")
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        print("****Alignment****")
        print(f"sequence: {alignment.title}")
        print(f"length: {alignment.length}")
        print(f"e-value: {hsp.expect}")
        print(f"score: {hsp.score}")
