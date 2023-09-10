# import fastaExt

push!(ARGS, "/lisc/user/pollak/software/fastaExt/tmp.faa")

build_header_index()
build_index()

append!(ARGS,
    ["PID-1232-NL-102.fixed_spades-metabat_cov.51.fa#NODE_507_length_42300_cov_2.150170_10#7583#8203#-1#00",
    "PID-1232-NL-102.fixed_spades-metabat_cov.48.fa#NODE_97_length_118217_cov_4.904535_89#107779#108495#-1#00",
    "PID-1232-NL-73.fixed_spades-metabat_cov.29.fa#NODE_547_length_38807_cov_3.559350_4#3650#4213#1#00"]
)

fastaext()

[pop!(ARGS) for _ in 1:3]
append!(ARGS,
    ["PID-1232-NL-102.fixed_spades-metabat_cov.51.fa", 
    "PID-1232-NL-73.fixed_spades-metabat_cov.29.fa"])

extgen()