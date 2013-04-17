#include <stdio.h>
#include "variant.h"
#include "cppmemory.h"
#include <Rinternals.h>


Variant::Variant(vector<Exon*>* exons)

{

	this->id = -1;

	this->name= "";

	//this->strand = strand;



	this->exonCount = exons->size();

	this->exons = new Exon*[exonCount];

	this->positions = new int[exonCount + 1];

	this->positions[0] = 1;



	/*this->codelen = (exonCount + 31) / 32;

	this->codes = new int[codelen];

	for (int i = 0; i < codelen; i++)

	{

	this->codes[i] = 0;

	}*/



	int i = 0;

	vector<Exon*>::const_iterator ei;

	for (ei = exons->begin(); ei != exons->end(); ei++)

	{

		Exon* exon = *ei;



		std::ostringstream out;

		out << (exon->id);

		(*this).name += out.str();

		(*this).name += ",";



		this->idmap[exon->id] = i;

		this->exons[i] = exon;

		this->positions[i + 1] = this->positions[i] + exon->length;

		//--- Variants have by default same strand as island. This must be changed for variants in islands with mixed strands. 
		//--- The case for known variants is solved in the importDataFrame function
		this->antisense = FALSE;

		/*int j = exon->num;

		this->codes[j / 32] |= 1 << (j % 32);*/

		i++;

	}



	(*this).name = (*this).name.substr(0, (*this).name.size()-1); //remove last ","

	this->length = this->positions[i] - 1;



	this->hashcode = gethash();

}



Variant::~Variant () {

	zaparray(exons); //delete [] exons;  //we delete the vector with exon pointers, not the exons themselves. After deleting a variant we want to keep the exons for other future variants

	zaparray(positions); //delete [] positions;

	//zaparray(codes); //delete [] codes;

}



int Variant::indexOf(int exonid)

{

	return this->idmap[exonid];

}



bool Variant::contains(Exon* e)

{

	return this->idmap.count(e->id) > 0;

}



bool Variant::contains(Fragment* frag)

{

  if (idmap.count(frag->left[0]) == 0) return false;

  for (int l = 1; l < frag->leftc; l++) {

    if ((idmap.count(frag->left[l]) == 0) || (this->indexOf(frag->left[l]) != 1+this->indexOf(frag->left[l-1]))) return false;

  }

  if (idmap.count(frag->right[0]) == 0) return false;

  for (int r = 1; r < frag->rightc; r++) {

    if ((idmap.count(frag->right[r]) == 0) || (this->indexOf(frag->right[r]) != 1+this->indexOf(frag->right[r-1]))) return false;

  }

  return true;

}



void Variant::toString(char *str)

{

	str[0] = '\0';



	for (int e = 0; e < exonCount; e++) sprintf(str, "%s,%i", str, exons[e]->id);



}



int Variant::compare(const Variant* other)

{

	if (this->exonCount < other->exonCount) 

	{

		return -1;

	} 

	else if (this->exonCount > other->exonCount) 

	{

		return +1;

	}



	for (int i = 0; i < exonCount; i++)

	{

		if (this->exons[i] != other->exons[i])

		{

			if (this->exons[i] > other->exons[i])

			{

				return +1;

			}

			// then (this->exons[i] < other->exons[i])

			return -1;

		}

	}

	/*

	for (int c = 0; c < this->codelen; c++) 

	{

		if (this->codes[c] < other->codes[c])

		{

            return -1;

		} 

		else if (this->codes[c] > other->codes[c]) 

		{

			return +1;

	    }

	}*/



	return 0;

}

int Variant::gethash()

{

	int h = 0;



	for (int i = 0; i < exonCount; i++)

	{

		Exon* exon = exons[i];

		h = exon->num + 17 * h;

		// hashing using a prime number

	}



	return h;



	/*int h = 0;



	for (int c = 0; c < codelen; c++)

	{

	h ^= codes[c];

	}



	return h;*/

}


