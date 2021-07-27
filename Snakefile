rule all: 
    input: 
        "snake_out_test.txt"
    
        
rule hello:
    output:
        "snake_out_test.txt"
    message: 
        "create the output file using a basic notebook"
    conda: 
        "conda_mutability_env.yml"
    shell: 
        "python snakefile_test.py"

