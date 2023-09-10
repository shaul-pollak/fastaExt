# fastaExt

A tool to index and quickly extract records from large fasta files written in the Julia programming language.

to compile the program run:
```bash
julia -t auto --project=. make.jl
cd compiled/lib  
ln -s julia/* .
```

for some reason PackageCompiler does not place the libraries in the correct dir, and this is why it is required to link them back to the lib directory.

## usage
To build an index for a fasta file in the current dir called prots.faa use
```bash
build_index prots.faa
```
this will create two files: prots.faa.idx1 and prots.faa.idx2 needed for pseudo random access into the large fasta file.

to print indivudual records from prots.faa to stdout after the index was created use
```bash
fastaExt prots.faa rec1 rec2 rec3 ...
```
where rec1, rec2, etc. are headers (what appears after the ">" sign in the fasta file)

if your fasta headers are formatted as ">genome#otherstuff" where a "#" character separates the genome name from the rest of the header, extracting all of the records belonging to a single genome can be achived by using
```bash
extgen prots.faa gen1 gen2 ...
```