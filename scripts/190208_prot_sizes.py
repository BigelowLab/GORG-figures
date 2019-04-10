from viruscope import readfa

filename = "/mnt/scgc/simon/simonsproject/gorg-tara_clustering/analyses/gorg_tara_tropics_80minid_m80_sizes.csv"
fasta = "/mnt/scgc/simon/simonsproject/gorg-tara_clustering/analyses/gorg_tara_tropics_80minid_m80.faa"

with open(filename, "w") as oh:
    for name, seq in readfa(open(fasta)):
        print(name, len(seq), file=oh, sep=",")
