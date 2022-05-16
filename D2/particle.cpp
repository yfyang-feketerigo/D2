#include "particle.h"
using std::iostream;
using std::string;
using std::vector;
Particle& seek_id(vector<Particle>& vec_p, size_t _id)
{

	for (size_t i = 0; i < vec_p.size(); i++)
	{
		if (vec_p[i].id == _id)
		{
			return vec_p[i];
		}
	}
	std::cerr << "particle " << _id << "not found!" << std::endl;
	throw std::runtime_error(("particle " + std::to_string(_id) + " not found!").c_str());

}

