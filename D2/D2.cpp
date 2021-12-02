#include "D2.h"
namespace D2
{
	void Configuration_neighbours::_update_neighbours()
	{
		const auto& vec_pa = this->get_particle();
		vec_neighbours.resize(vec_pa.size());
		for (size_t i = 0; i < vec_pa.size(); i++)
		{
			vec_neighbours[i].p_center_pa = &vec_pa[i];
			vec_neighbours[i].pvec_neighbours.clear();
		}
		for (size_t i = 0; i < vec_pa.size(); i++)
		{
			const Particle& pai = vec_pa[i];
			for (size_t j = i + 1; j < vec_pa.size(); j++)
			{
				const Particle& paj = vec_pa[j];
				VectorDd rij_PBC = get_PBC_rij(pai, paj);
				bool flag_compute_drij = true;
				for (size_t k = 0; k < DIMENSION; k++) // pre-check distance
				{
					flag_compute_drij = (flag_compute_drij && (abs(rij_PBC[k]) < r_cut));
				}
				if (flag_compute_drij)
				{
					double drij = sqrt(rij_PBC.dot(rij_PBC));
					if (drij < r_cut)
					{
						vec_neighbours[i].pvec_neighbours.push_back(&paj);
						vec_neighbours[j].pvec_neighbours.push_back(&pai);
					}
				}
			}
		}
		flag_neighbous_update = true;
	}

	const Neighbours& Configuration_neighbours::get_neighbours(size_t pid) const
	{
		for (size_t i = 0; i < vec_neighbours.size(); i++)
		{
			if (pid == vec_neighbours[i].p_center_pa->id)
			{
				return vec_neighbours[i];
			}
		}
		std::string err_message = "particle id " + std::to_string(pid) + " not found!\n";
		throw(std::exception(err_message.c_str()));
		return vec_neighbours.front();
	}

	VectorDd Configuration_neighbours::to_box_coordination(const Particle& pa) const
	{
		VectorDd pa_vector_cartesian;
		pa_vector_cartesian(0) = pa.rx;
		pa_vector_cartesian(1) = pa.ry;
		pa_vector_cartesian(2) = pa.rz;
		auto pa_vector_box = m_vector_cartesian_to_box * pa_vector_cartesian;
		return pa_vector_box;
	}

	VectorDd Configuration_neighbours::get_PBC_rij(const Particle& pai, const Particle& paj) const
	{
		VectorDd pai_vector_box = to_box_coordination(pai);
		VectorDd paj_vector_box = to_box_coordination(paj);
		VectorDd rij_box = pai_vector_box - paj_vector_box;

		for (size_t i = 0; i < rij_box.size(); i++)
		{
			rij_box[i] -= round(rij_box[i]);
		}
		;
		VectorDd rij_cartesian = m_vector_box_to_cartesian * rij_box;
		return rij_cartesian;
	}


	void D2::compute_Xij()
	{
		Xij.setZero();
		const Particle& pa_0_t = *p_pa;
		auto& pa_0_t_neighbours = config_t->get_neighbours(pa_0_t);
		const Particle& pa_0_t_ = config_t_->get_particle(pa_0_t.id);
		for (size_t n = 0; n < pa_0_t_neighbours.pvec_neighbours.size(); n++)
		{
			const Particle& pa_n_t = *pa_0_t_neighbours.pvec_neighbours[n];
			const Particle& pa_n_t_ = config_t_->get_particle(pa_n_t.id);
			VectorDd rn0_t = config_t->get_PBC_rij(pa_n_t, pa_0_t);
			VectorDd rn0_t_ = config_t_->get_PBC_rij(pa_n_t_, pa_0_t_);
			Xij += rn0_t * rn0_t_.transpose();
		}
		flag_computed_Xij = true;
		return;
	}

	void D2::compute_Yij()
	{
		Yij.setZero();
		const Particle& pa_0_t = *p_pa;
		auto& pa_0_t_neighbours = config_t->get_neighbours(pa_0_t);
		const Particle& pa_0_t_ = config_t_->get_particle(pa_0_t.id);
		for (size_t n = 0; n < pa_0_t_neighbours.pvec_neighbours.size(); n++)
		{
			const Particle& pa_n_t = *pa_0_t_neighbours.pvec_neighbours[n];
			const Particle& pa_n_t_ = config_t_->get_particle(pa_n_t.id);
			VectorDd rn0_t_ = config_t_->get_PBC_rij(pa_n_t_, pa_0_t_);
			Yij += rn0_t_ * rn0_t_.transpose();
		}
		flag_computed_Yij = true;
		return;
	}


	void D2::compute_epsilonij()
	{
		epsilonij.setZero();
		compute_Xij();
		compute_Yij();
		epsilonij = Xij * (Yij.inverse().transpose()) - DELTAij; //DELTAij = 1,0,0
																 //          0,1,0
																 //          0,0,1
		flag_computed_epsilonij = true;
		return;
	}

	void D2::compute_D2()
	{
		compute_epsilonij();
		const Particle& pa_0_t = *p_pa;
		const Particle& pa_0_t_ = config_t_->get_particle(pa_0_t.id);
		const auto& pvec_ngbrs_0_t = config_t->get_neighbours(pa_0_t).pvec_neighbours;
		d2 = 0;
		for (size_t n = 0; n < pvec_ngbrs_0_t.size(); n++)
		{
			const Particle& pa_n_t = *pvec_ngbrs_0_t[n];
			const Particle& pa_n_t_ = config_t_->get_particle(pa_n_t.id);
			VectorDd rn0_t = config_t->get_PBC_rij(pa_n_t, pa_0_t);
			VectorDd rn0_t_ = config_t_->get_PBC_rij(pa_n_t_, pa_0_t_);
			VectorDd Dn_vector = rn0_t - (DELTAij + epsilonij) * rn0_t_;
			double Dn = Dn_vector.dot(Dn_vector);
			d2 += Dn;
		}
		d2 /= pvec_ngbrs_0_t.size();
		flag_computed_d2 = true;
	}
}
