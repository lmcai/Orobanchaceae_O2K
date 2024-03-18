import tkinter as tk
from tkinter import filedialog
from Bio import SeqIO

class SequenceRetrieverApp:
    def __init__(self, master):
        self.master = master
        self.master.title("Sequence Retriever")

        self.fasta_label = tk.Label(master, text="Select FASTA file:")
        self.fasta_label.grid(row=0, column=0, sticky="w", padx=10, pady=5)

        self.fasta_entry = tk.Entry(master, width=50)
        self.fasta_entry.grid(row=0, column=1, padx=10, pady=5)

        self.fasta_button = tk.Button(master, text="Browse", command=self.browse_fasta)
        self.fasta_button.grid(row=0, column=2, padx=5, pady=5)

        self.ids_label = tk.Label(master, text="Enter IDs (separated by commas):")
        self.ids_label.grid(row=1, column=0, sticky="w", padx=10, pady=5)

        self.ids_entry = tk.Entry(master, width=50)
        self.ids_entry.grid(row=1, column=1, padx=10, pady=5)
        self.ids_entry.bind('<FocusOut>', self.retrieve_sequences)

        self.retrieve_button = tk.Button(master, text="Retrieve Sequences", command=self.retrieve_sequences)
        self.retrieve_button.grid(row=2, column=1, pady=10)

        self.output_text = tk.Text(master, height=10, width=50)
        self.output_text.grid(row=3, columnspan=3, padx=10, pady=5)

    def browse_fasta(self):
        fasta_file = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta"), ("All files", "*.*")])
        self.fasta_entry.delete(0, tk.END)
        self.fasta_entry.insert(0, fasta_file)

    def retrieve_sequences(self, event=None):
        fasta_file = self.fasta_entry.get()
        ids = [x.strip() for x in self.ids_entry.get().split(',') if x.strip()]
        
        sequences = {}
        with open(fasta_file, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences[record.id] = record.seq

        output = ""
        for id in ids:
            if id in sequences:
                output += f">{id}\n{sequences[id]}\n"
            else:
                output += f"Sequence with ID '{id}' not found!\n"

        self.output_text.delete(1.0, tk.END)
        self.output_text.insert(tk.END, output)

def main():
    root = tk.Tk()
    app = SequenceRetrieverApp(root)
    root.mainloop()

if __name__ == "__main__":
    main()
