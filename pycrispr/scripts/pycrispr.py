#!/usr/bin/env python3
import click

#command line parser
@click.command()
@click.option("-l","--library", required=True, help="CRISPR-Cas9 library")
@click.option("-m","--mismatch", default=0, show_default=True, type=int, help="Number of mismatches allowed during alignment")

def cli(library,mismatch):
	 """CRISPR-Cas9 screen analysis"""
	 click.echo(f"{library} library selected")
	 click.echo(f"Mismatches allowed: {mismatch}")




