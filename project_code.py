import tkinter as tk
from tkinter import messagebox
from Bio import Entrez  # type: ignore
from Bio import SeqIO  # type: ignore
from Bio.Blast import NCBIWWW #type: ignore
from Bio.Blast import NCBIXML #type: ignore

Entrez.email = "youremail@abc.com"

def fet_data(event=None):
    
    try:
        acc_number = ent.get().strip()
        handle = Entrez.efetch(db="nucleotide", id=acc_number, rettype="fasta", retmode="text")
        record = handle.read()
        text_box.delete("1.0", tk.END)   
        text_box.insert(tk.END, record) 
    except Exception:
        messagebox.showinfo("info", "Enter valid accession ID!")

def gc_count():
    
    try:
        acc_number = ent.get().strip()
        handle = Entrez.efetch(db="nucleotide", id=acc_number, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        sequence = record.seq
        GC = 100 * float(sequence.count("G") + sequence.count("C")) / len(sequence)
        text_box.delete("1.0", tk.END)   
        text_box.insert(tk.END, f"GC Content: {GC:.2f}%") 
    except Exception:
        messagebox.showinfo("info", "Enter valid accession ID!")

def c_d():
    
    try:
        acc_number = ent.get().strip()
        handle = Entrez.efetch(db="nucleotide", id=acc_number, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        sequence = record.seq
        Transcribe= sequence.transcribe()
        Translate= sequence.translate()
        text_box.delete("1.0", tk.END)   
        text_box.insert(tk.END,f"DNA: {sequence}\n\n RNA: {Transcribe}\n\nPROTEIN: {Translate}") 
    except Exception:
        messagebox.showinfo("info", "Enter valid accession ID!")

def blast():
    try:
        acc_number = ent.get().strip()

       
        handle = Entrez.efetch(db="nucleotide", id=acc_number,
                               rettype="fasta", retmode="text")
        record = handle.read()
        nucleotide_seq = record

       
        result_handle_n = NCBIWWW.qblast("blastn", "nt", nucleotide_seq)
        blast_records = NCBIXML.parse(result_handle_n)

        from textwrap import wrap  

        output_text = f"BLAST results for {acc_number}:\n\n"

        for blast_record in blast_records:
            count = 0
            for alignment in blast_record.alignments:
                if count >= 3:  
                    break
                output_text += f"Sequence: {alignment.title}\n"
                output_text += f"Length: {alignment.length}\n"

                for hsp in alignment.hsps:
                    output_text += f"Score: {hsp.score}, E-value: {hsp.expect}\n"

                    query_lines = wrap(hsp.query, 60)
                    match_lines = wrap(hsp.match, 60)
                    sbjct_lines = wrap(hsp.sbjct, 60)

                    for q, m, s in zip(query_lines, match_lines, sbjct_lines):
                        output_text += f"Query:   {q}\n"
                        output_text += f"         {m}\n"
                        output_text += f"Subject: {s}\n\n"

                    break  
                count += 1
            break  

        text_box.delete("1.0", tk.END)
        text_box.insert(tk.END, output_text)

    except Exception as e:
        messagebox.showinfo("Error", f"Enter valid accession ID!\n{e}")

root = tk.Tk()
root.title("GeneScope")
root.geometry("1000x600+50+50")


lab1=tk.Label(root,text="Accession id")
lab1.pack()
ent = tk.Entry(root, width=40)
ent.pack(pady=20, padx=20)
ent.focus_set()

frame = tk.Frame(root,)
frame.pack(fill="both", expand=True, padx=10, pady=10)


scrollbar = tk.Scrollbar(frame)
scrollbar.pack(side="right", fill="y")

text_box = tk.Text(frame, wrap="word", yscrollcommand=scrollbar.set)
text_box.pack(side="left", fill="both", expand=True)

scrollbar.config(command=text_box.yview)

button_frame = tk.Frame(root)
button_frame.pack(pady=20)

but = tk.Button(button_frame, text="Fetch Sequence", command=fet_data, width=20)
but.grid(row=0, column=0, padx=10, pady=5)

but1 = tk.Button(button_frame, text="GC Content", command=gc_count, width=20)
but1.grid(row=0, column=1, padx=10, pady=5)

but2 = tk.Button(button_frame, text="Central Dogma", command=c_d, width=20)
but2.grid(row=1, column=0, padx=10, pady=5)

but3 = tk.Button(button_frame, text="Blast", command=blast, width=20)
but3.grid(row=1, column=1, padx=10, pady=5)

root.bind("<Return>", fet_data)

root.mainloop()

