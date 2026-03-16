import shlex
import subprocess
import os
from string import strip


mainpath = "/home/BMMM_training/defazio/tirocinio_postlaurea/itsonedb"

#TEST SCRIPT HMMSEARCH
program = os.path.join(mainpath, "scripts/paral_estrazione_localizzazione_v1.py")
folder = os.path.join(mainpath, "hmmsearch_prova/debug")
out = os.path.join(mainpath, "hmmsearch_prova/debug/hmmsearch_test.csv")
log = os.path.join(mainpath, "hmmsearch_prova/debug/hmmsearch_test.log")
command = "python2 %s -f %s -o %s" % (program, folder, out)
with open(log, "w") as log_out:
    subprocess.call(shlex.split(command), stdout=log_out)

hmmsearch_test = {
    "prova_r131\tBZ470106\tN\t2.3e-11\t0.03\t0.96\t1\t49\t1.9e-25\t0.61\t0.98\t315\t409\n": """TEST1
    UN GENE 18S E UN GENE 5.8S STRAND FORWARD ------- > OK""",
    "prova_r131\tCT573131\tN\t0.0\t1.0\t0.88\t63677\t65612\t8.7e-27\t0.99\t0.92\t66007\t66161\n" : """TEST2
    UN GENE 18S CON MARGINE DESTRO VALIDATO, IL SUCCESSIVO CON MARGINE SINISTRO VALIDATO E 
    IL GENE 5.8S TRA I DUE --------- > OK""",
    "prova_r131\tTEST3\tN\t0.0\t1.0\t0.88\t303\t2103\t8.7e-27\t0.99\t0.92\t2253\t2403\n" : """TEST3
    DUE GENI 18S CON MARGINE DESTRO VALIDATO, 5.8S DEL 1* GENE ------- > OK""",
    "prova_r131\tTEST4\tN\t0.0\t1.0\t0.88\t2553\t4353\t8.7e-27\t0.99\t0.92\t4503\t4656\n": """TEST4
    DUE GENI 18S CON MARGINE DESTRO VALIDATO, 5.8S DEL 2* GENE ------- > OK"""
}

with open(log, "r") as log_out:
    less_log = log_out.readlines()
    for line in less_log:
        if "Traceback" in line:
            print "HMMSEARCH SCRIPT ERROR"

with open(out, "r") as out_out:
    less_out = out_out.readlines()
    print "TEST SCRIPT HMMSEARCH"
    for line in less_out[1:]:
        if line in hmmsearch_test.keys():
            print hmmsearch_test[line]
        else:
            s = map(strip, line.split("\t"))
            print "%s FALLITO" % s[1]

os.remove(log)
os.remove(out)

#TEST SCRIPT NHMMER
program = os.path.join(mainpath, "scripts/paral_estrazione_localizzazione_v2.py")
folder = os.path.join(mainpath, "nhmmer_prova/debug")
out = os.path.join(mainpath, "nhmmer_prova/debug/nhmmer_test.csv")
log = os.path.join(mainpath, "nhmmer_prova/debug/nhmmer_test.log")
command = "python2 %s -f %s -o %s" % (program, folder, out)
with open(log, "w") as log_out:
    subprocess.call(shlex.split(command), stdout=log_out)

nhmmer_tests = {
    "TEST1": ("prova_r131\tTEST1\t410\t4.6e-08\t0.03\t0.96\t1\t49\t5.8e-22\t0.61\t0.98\t315\t409", """TEST1
     UN GENE 18S E UN GENE 5.8S STRAND FORWARD -------- > OK"""),
    "TEST2": ("prova_r131\tTEST2\t762\t4.9e-35\t0.09\t0.88\t762\t590\t1.2e-32\t0.98\t0.93\t324\t171", """TEST2
     UN GENE 18S E UN GENE 5.8S STRAND REVERSE -------- > OK"""),
    "TEST3": ("prova_r131\tTEST3\t76958\t0.0\t1.0\t0.87\t63677\t65613\t1.5e-23\t0.99\t0.92\t66007\t66161", """TEST3
     UN GENE 18S CON MARGINE DESTRO VALIDATO, IL SUCCESSIVO CON MARGINE SINISTRO VALIDATO E 
     IL GENE 5.8S TRA I DUE. STRAND FORWARD ----------- > OK"""),
    "TEST4": ("prova_r131\tTEST4\t136866\t0.0\t1.0\t0.87\t13300\t11364\t1.5e-23\t0.99\t0.92\t10970\t10816","""TEST4 
     UN GENE 18S CON MARGINE DESTRO VALIDATO, IL SUCCESSIVO CON MARGINE SINISTRO VALIDATO E 
     IL GENE 5.8S TRA I DUE. STRAND REVERSE ----------- > OK"""),
    "TEST5": ("prova_r131\tTEST5\t485\t1.1e-09\t0.03\t0.95\t51\t100\t9.3e-27\t0.75\t0.94\t301\t416","""TEST5 
     DUE GENI 18S CON MARGINE DESTRO VALIDATO, 5.8S DEL 2* GENE. STRAND FORWARD ----- > OK"""),
    "TEST6": ("prova_r131\tTEST6\t777\t1.1\t0.03\t0.88\t109\t161\t6.2e-22\t0.82\t0.92\t398\t523", """TEST6 
     DUE GENI 18S CON MARGINE DESTRO VALIDATO, 5.8S DEL 1* GENE. STRAND FORWARD ----- > OK"""),
    "TEST7": ("prova_r131\tTEST7\t777\t2.7e-23\t0.01\t0.85\t4352\t2553\t6.2e-22\t0.99\t0.92\t2403\t2253","""TEST7
     DUE GENI 18S CON MARGINE DESTRO VALIDATO, 5.8S DEL 1* GENE. STRAND REVERSE ----- > OK"""),
    "TEST8": ("prova_r131\tTEST8\t777\t2.7e-23\t0.01\t0.85\t2103\t303\t6.2e-22\t0.99\t0.92\t153\t1", """TEST8
     DUE GENI 18S CON MARGINE DESTRO VALIDATO, 5.8S DEL 2* GENE. STRAND REVERSE ----- > OK""")
}

with open(log, "r") as log_out:
    less_log = log_out.readlines()
    for line in less_log:
        if "Traceback" in line:
            print "NHMMER SCRIPT ERROR"

with open(out, "r") as out_out:
    less_out = out_out.readlines()
    print "TEST SCRIPT NHMMER"
    for line in less_out[1:]:
        s = line.split("\t")
        p = nhmmer_tests[s[1]][0].split("\t")
        if (s[6] == p[6]) and (s[7] == p[7]) and (s[11] == p[11]) and (s[12] == p[12]):
            print nhmmer_tests[s[1]][1]
        else:
            print "%s FALLITO" % s[1]
#os.remove(log)
#os.remove(out)
