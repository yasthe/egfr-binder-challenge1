#!/usr/bin/env python3
"""Bringt submission/ auf den Stand dessen, was tatsaechlich eingereicht wurde.

submission/designs.csv enthaelt danach nur die eingereichten Designs;
submission/metrics.csv behaelt alle sieben und bekommt eine Spalte "submitted",
damit nachvollziehbar bleibt, welche gepruefte Kandidaten es nicht wurden.

Aufruf aus dem Repositorium heraus:  python scripts/mark_submitted.py
"""
import os, sys
import pandas as pd

EINGEREICHT = {"EGFRpH_C_l55_s661569_mpnn4",
               "EGFRpH_C_l55_s661569_mpnn6",
               "EGFRpH_C_l45_s295112_mpnn1"}
GRUND = "novelty score 2/4, below the platform threshold of 3/4"

d = pd.read_csv("submission/designs.csv")
raus = sorted(set(d["name"]) - EINGEREICHT)
fehlt = sorted(EINGEREICHT - set(d["name"]))
if fehlt:
    sys.exit(f"Nicht in designs.csv gefunden: {fehlt}")

d[d["name"].isin(EINGEREICHT)].to_csv("submission/designs.csv", index=False)
with open("submission/designs.fasta") as fh:
    bloecke = fh.read().split(">")
with open("submission/designs.fasta", "w") as fh:
    for b in bloecke:
        if b.strip() and b.split("\n")[0].strip() in EINGEREICHT:
            fh.write(">" + b)

m = pd.read_csv("submission/metrics.csv")
sp = "Design" if "Design" in m.columns else m.columns[0]
m["submitted"] = m[sp].isin(EINGEREICHT)
m["withheld_reason"] = [""if x else GRUND for x in m["submitted"]]
m.to_csv("submission/metrics.csv", index=False)

print(f"designs.csv: {len(EINGEREICHT)} eingereichte Designs")
print(f"zurueckgehalten: {', '.join(raus) if raus else 'keine'}")
print("metrics.csv: Spalten submitted und withheld_reason ergaenzt")
