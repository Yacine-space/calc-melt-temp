from flask import Flask, request, render_template
from Bio.Seq import Seq
import math

app = Flask(__name__)


def analyse_sequence(sequence, conc_na=50.0):
    sequence = sequence.upper().strip()
    sequence_list = []
    for n in sequence:
        if n == "A" or n == "C" or n == "G" or n == "T":
            sequence_list.append(n)
        else:
            return None

    s1 = Seq(sequence)
    complement = s1.reverse_complement()
    num_a = sequence_list.count("A")
    num_t = sequence_list.count("T")
    num_c = sequence_list.count("C")
    num_g = sequence_list.count("G")
    longueur = num_c + num_g + num_a + num_t
    pourc_gc = (num_c + num_g) / longueur * 100
    pourc_gc = round(pourc_gc, 1)
    nb_gc = num_g + num_c
    nb_at = num_a + num_t

    tm_val_basic = 64.9 + 41 * ((nb_gc-16.4)/longueur)
    tm_val_salt_adjusted = 100.5 + (41 * (nb_gc / longueur)) - 820 / longueur + 16.6 * math.log10(conc_na / 1000)
    tm_val_basic = round(tm_val_basic,1)
    tm_val_salt_adjusted = round(tm_val_salt_adjusted,1)
    hyb_temp = tm_val_salt_adjusted-5
    sequence_colored = sequence_coloree(sequence)
    complement_colored = sequence_coloree(str(complement))
    return {
        "longueur_seq": longueur,
        "nombre_A": num_a,
        "nombre_C": num_c,
        "nombre_G": num_g,
        "nombre_T": num_t,
        "pourcentage_GC": pourc_gc,
        "Tm_basic": tm_val_basic,
        "Tm_adjusted": tm_val_salt_adjusted,
        "sequence": sequence,
        "complement": complement,
        "sequence_colored": sequence_colored,
        "complement_colored": complement_colored,
        "conc_na": conc_na,
        "temperature_hybridation": hyb_temp
    }


def read_fasta(file):
    sequence = ""
    for line in file:
        decoded_line = line.decode("utf-8").strip()
        if not decoded_line.startswith(">"):
            sequence += decoded_line
    return sequence


def sequence_coloree(seq):
    colored = ""
    for base in seq.upper():
        colored += f'<span class="base {base}">{base}</span>'
    return colored

@app.route("/", methods=["GET", "POST"])
def home():
    resultat = None

    if request.method == "POST":
        conc_na_raw = request.form.get("conc_na", "50").strip()
        try:
            conc_na = float(conc_na_raw)
            if conc_na <= 0:
                conc_na = 50.0
        except ValueError:
            conc_na = 50.0

        sequence = ""
        if request.form.get("sequence", "").strip():
            sequence = request.form["sequence"]
        elif "file" in request.files and request.files["file"].filename != "":
            sequence = read_fasta(request.files["file"])

        if sequence:
            resultat = analyse_sequence(sequence, conc_na)

    return render_template("index.html", resultat=resultat)

if __name__ == "__main__":
    app.run(debug=True)
