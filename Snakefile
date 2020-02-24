#!python


# rev: revisions
include: "snakefiles/rev/Snakefile_read_sampling"
include: "snakefiles/rev/Snakefile_read_sampling_pmd"
include: "snakefiles/rev/Snakefile_pca"
include: "snakefiles/rev/Snakefile_eigenstrat"
include: "snakefiles/rev/Snakefile_admixture"
include: "snakefiles/rev/Snakefile_y"
include: "snakefiles/rev/Snakefile_damage"
include: "snakefiles/rev/Snakefile_dystruct"
include: "snakefiles/ExtractEigenstrat/Snakefile"
include: "snakefiles/rev/Snakefile_read_sampling_xy"


rule none:
    input: 'Snakefile'
    run: print("ancient-sardinia")
