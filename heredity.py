import pdb
import csv
import itertools
import sys

PROBS = {

	# Unconditional probabilities for having gene
	"gene": {
		2: 0.01,
		1: 0.03,
		0: 0.96
	},

	"trait": {

		# Probability of trait given two copies of gene
		2: {
			True: 0.65,
			False: 0.35
		},

		# Probability of trait given one copy of gene
		1: {
			True: 0.56,
			False: 0.44
		},

		# Probability of trait given no gene
		0: {
			True: 0.01,
			False: 0.99
		}
	},

	# Mutation probability
	"mutation": 0.01
}


def main():

	# Check for proper usage
	if len(sys.argv) != 2:
		sys.exit("Usage: python heredity.py data.csv")
	people = load_data(sys.argv[1])

	# Keep track of gene and trait probabilities for each person
	probabilities = {
		person: {
			"gene": {
				2: 0,
				1: 0,
				0: 0
			},
			"trait": {
				True: 0,
				False: 0
			}
		}
		for person in people
	}

	# Loop over all sets of people who might have the trait
	names = set(people)
	for have_trait in powerset(names):

		# Check if current set of people violates known information
		fails_evidence = any(
			(people[person]["trait"] is not None and
			 people[person]["trait"] != (person in have_trait))
			for person in names
		)
		if fails_evidence:
			continue

		# Loop over all sets of people who might have the gene
		for one_gene in powerset(names):
			for two_genes in powerset(names - one_gene):
				# Update probabilities with new joint probability
				p = joint_probability(people, one_gene, two_genes, have_trait)
				update(probabilities, one_gene, two_genes, have_trait, p)

	# Ensure probabilities sum to 1
	normalize(probabilities)

	# Print results
	for person in people:
		print(f"{person}:")
		for field in probabilities[person]:
			print(f"  {field.capitalize()}:")
			for value in probabilities[person][field]:
				p = probabilities[person][field][value]
				print(f"    {value}: {p:.4f}")


def load_data(filename):
	"""
	Load gene and trait data from a file into a dictionary.
	File assumed to be a CSV containing fields name, mother, father, trait.
	mother, father must both be blank, or both be valid names in the CSV.
	trait should be 0 or 1 if trait is known, blank otherwise.
	"""
	data = dict()
	with open(filename) as f:
		reader = csv.DictReader(f)
		for row in reader:
			name = row["name"]
			data[name] = {
				"name": name,
				"mother": row["mother"] or None,
				"father": row["father"] or None,
				"trait": (True if row["trait"] == "1" else
						  False if row["trait"] == "0" else None)
			}
	return data


def powerset(s):
	"""
	Return a list of all possible subsets of set s.
	"""
	s = list(s)
	return [
		set(s) for s in itertools.chain.from_iterable(
			itertools.combinations(s, r) for r in range(len(s) + 1)
		)
	]


def joint_probability(people, one_gene, two_genes, have_trait):
	"""
	Compute and return a joint probability.

	The probability returned should be the probability that
			* everyone in set `one_gene` has one copy of the gene, and
			* everyone in set `two_genes` has two copies of the gene, and
			* everyone not in `one_gene` or `two_gene` does not have the gene, and
			* everyone in set `have_trait` has the trait, and
			* everyone not in set` have_trait` does not have the trait.
	"""

	final_prob = 1

	for person in people:
		probability = 1
        
		# probability that the person has one gene
		if person in one_gene:
			# find their parents
			mom = people[person]['mother']
			dad = people[person]['father']
			# if person has no parents
			if mom == None and dad == None:
				# unconditional probability of having one gene
				probability = PROBS["gene"][1]
			# if person has parents
			else:
				# gene comes from mother and not father (probability)
				if mom in one_gene:
					# mom has one harmful gene
					# harmful gene is selected from mom with p = 0.5
					probability = 0.5
				elif mom in two_genes:
					# mom has both harmful genes
					probability = 1 - PROBS["mutation"]
				else:
					# mom has no harmful genes
					probability = PROBS["mutation"]
				
				if dad in one_gene:
					# harmless gene is selected from dad with p = 0.5
					probability = (probability * 0.5) + ((1 - probability) * 0.5)
					# (multiplying because joint probability)
				elif dad in two_genes:
					# no matter which gene is selected
					# it mutates into the harmless form with p = PROBS["mutation"]
					probability = (probability * PROBS["mutation"]) + ((1 - PROBS["mutation"]) * (1 -probability))
				else:
					# this means dad has no harmful gene
					# no matter which harmless gene is selected
					# it stays harmless with probability 1 - PROBS["mutation"]
					probability = (probability * (1- PROBS["mutation"])) + ((1 -probability) * PROBS["mutation"])



			# given one harmful gene
			# probability of trait | not
			probability *= PROBS["trait"][1][person in have_trait]

		elif person in two_genes:
			# find their parents
			mom = people[person]['mother']
			dad = people[person]['father']
			# if a person has no parents
			if mom == None and dad == None:
				# unconditional probability of having two genes
				probability = PROBS["gene"][2]
			else:
				# both genes harmful
				# one harmful gene comes from each parent
				if mom in one_gene:
					probability = 0.5
				elif mom in two_genes:
					probability = 1 - PROBS["mutation"]
				else:
					probability = PROBS["mutation"]

				if dad in one_gene:
					probability = (probability * 0.5) + ((1 - probability) * 0.5)
				elif dad in two_genes:
					probability = (probability * PROBS["mutation"]) + ((1 - PROBS["mutation"]) * (1 -probability))
				else:
					probability = (probability * (1- PROBS["mutation"])) + ((1 -probability) * PROBS["mutation"])

			# given two genes,
			# probablity of trait
			probability *= PROBS["trait"][2][person in have_trait]

		else:
			# person has no harmful genes
			# find his parents
			mom = people[person]['mother']
			dad = people[person]['father']

			if mom == None and dad == None:
				# unconditional probability of having no harmful genes
				probability = PROBS["gene"][0]
			else:
				# no harmful genes come from parents
				if mom in one_gene:
					probability = 0.5
				elif mom in two_genes:
					probability = 1 - PROBS["mutation"]
				else:
					probability = PROBS["mutation"]

				if dad in one_gene:
					probability = (probability * 0.5) + ((1 - probability) * 0.5)
				elif dad in two_genes:
					probability = (probability * PROBS["mutation"]) + ((1 - PROBS["mutation"]) * (1 -probability))
				else:
					probability = (probability * (1- PROBS["mutation"])) + ((1 -probability) * PROBS["mutation"])

			# given no harmful genes, probability of trait
			probability *= PROBS["trait"][0][person in have_trait]

		final_prob *= probability
	
	return final_prob

def update(probabilities, one_gene, two_genes, have_trait, p):
	"""
	Add to `probabilities` a new joint probability `p`.
	Each person should have their "gene" and "trait" distributions updated.
	Which value for each distribution is updated depends on whether
	the person is in `have_gene` and `have_trait`, respectively.
	"""
	
	for person in probabilities:
		if person in one_gene:
			probabilities[person]["gene"][1] += p	
		elif person in two_genes:
			probabilities[person]["gene"][2] += p
		else:
			probabilities[person]["gene"][0] += p

		probabilities[person]["trait"][person in have_trait] += p

def normalize(probabilities):
    """
    Update `probabilities` such that each probability distribution
    is normalized (i.e., sums to 1, with relative proportions the same).
    """
	

    for person in probabilities:
        total_genes = sum(probabilities[person]["gene"].values())
        for i in probabilities[person]["gene"]:
            probabilities[person]["gene"][i] = probabilities[person]["gene"][i] / total_genes

        total_traits = sum(probabilities[person]["trait"].values())
        for i in probabilities[person]["trait"]:
            probabilities[person]["trait"][i] = probabilities[person]["trait"][i] / total_traits
	

if __name__ == "__main__":
    main()
