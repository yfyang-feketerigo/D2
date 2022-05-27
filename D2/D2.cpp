#include "D2.h"
namespace D2
{
	void Configuration_neighbours::_update_neighbours()
	{
		vec_neighbours.clear();
		double rho = get_particle_num() / get_lx() / get_ly() / get_lz();
		size_t reserve_neighbour_size = size_t(rho * rcut * rcut * rcut) + NEIGHBOUR_SIZE_REDUNDANCY;

		const auto& vec_pa = this->get_particle();
		vec_neighbours.resize(vec_pa.size());
		for (size_t i = 0; i < vec_pa.size(); i++)
		{
			vec_neighbours[i].p_center_pa = &vec_pa[i];
			vec_neighbours[i].vec_neighbours.reserve(reserve_neighbour_size);
			vec_neighbours[i].vec_neighbours.clear();
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
					flag_compute_drij = (flag_compute_drij && (abs(rij_PBC[k]) < rcut));
				}
				if (flag_compute_drij)
				{
					double drij = sqrt(rij_PBC.dot(rij_PBC));
					if (drij < rcut)
					{
#ifndef NO_STORE_NEIGHBOUR_RIJ // store rij information
						Neighbour j_in_i{ &paj,rij_PBC,drij };
						Neighbour i_in_j{ &pai,-rij_PBC,drij };
#else // NOT store rij information
						Neighbour j_in_i{ &paj };
						Neighbour i_in_j{ &pai };
#endif // !NO_STORE_NEIGHBOUR_RIJ
						vec_neighbours[i].vec_neighbours.push_back(j_in_i);
						vec_neighbours[j].vec_neighbours.push_back(i_in_j);
					}
				}
			}
		}
		flag_neighbours_update = true;
		flag_neighbours_sorted = false;
	}

	Configuration_neighbours::Configuration_neighbours(const Configuration_neighbours& config) :Configuration::Configuration(config)
	{
		this->rcut = config.rcut;
		this->m_base_cartesian_to_box = config.m_base_cartesian_to_box;
		this->m_base_box_to_cartesian = config.m_base_box_to_cartesian;
		this->m_vector_cartesian_to_box = this->m_base_box_to_cartesian;
		this->m_vector_box_to_cartesian = this->m_base_cartesian_to_box;
		if (config.flag_neighbours_update)
		{
			this->_update_neighbours();
		}
		else
		{
			this->flag_neighbours_update = false;
			this->vec_neighbours = std::move(std::vector<Neighbours>());
		}
		if (config.flag_neighbours_sorted)
		{
			this->sort_neighbours_as_center_pid();
		}
		else
		{
			this->flag_neighbours_sorted = false;
			this->pvec_neighbours_sorted = std::move(std::vector<const Neighbours*>());
		}

	}

	Configuration_neighbours& Configuration_neighbours::operator=(const Configuration_neighbours& config)
	{
		if (this == &config)
		{
			return *this;
		}
		Configuration::Configuration:: operator= (config);
		this->rcut = config.rcut;
		this->m_base_cartesian_to_box = config.m_base_cartesian_to_box;
		this->m_base_box_to_cartesian = config.m_base_box_to_cartesian;
		this->m_vector_cartesian_to_box = this->m_base_box_to_cartesian;
		this->m_vector_box_to_cartesian = this->m_base_cartesian_to_box;
		if (config.flag_neighbours_update)
		{
			this->_update_neighbours();
		}
		else
		{
			this->flag_neighbours_update = false;
			this->vec_neighbours = std::move(std::vector<Neighbours>());
		}
		if (config.flag_neighbours_sorted)
		{
			this->sort_neighbours_as_center_pid();
		}
		else
		{
			this->flag_neighbours_sorted = false;
			this->pvec_neighbours_sorted = std::move(std::vector<const Neighbours*>());
		}

		return *this;
		// TODO: 在此处插入 return 语句
	}

	void Configuration_neighbours::sort_neighbours_as_center_pid()
	{
		pvec_neighbours_sorted.resize(vec_neighbours.size());
		for (size_t i = 0; i < pvec_neighbours_sorted.size(); i++)
		{
			pvec_neighbours_sorted[i] = &vec_neighbours[i];
		}
		auto comp_cid = [](const Neighbours* a, const Neighbours* b) {return a->p_center_pa->id < b->p_center_pa->id; };
		std::sort(pvec_neighbours_sorted.begin(), pvec_neighbours_sorted.end(), comp_cid);
		flag_neighbours_sorted = true;
	}

	const Neighbours& Configuration_neighbours::get_neighbours(size_t pid) const
	{
		if (flag_neighbours_sorted == false)
		{
			for (size_t i = 0; i < vec_neighbours.size(); i++)
			{
				if (pid == vec_neighbours[i].p_center_pa->id)
				{
					return vec_neighbours[i];
				}
			}
			std::string err_message = "particle id " + std::to_string(pid) + " not found!\n";
			throw(std::runtime_error(err_message.c_str()));
			return vec_neighbours.front();
		}
		else
		{
#ifdef SAFE_GET_NEIGHBOURS
			if (pid <= 0 || pid > pvec_neighbours_sorted.size())
			{
				throw std::runtime_error(("illegal pid " + std::to_string(pid) + " when get_neighbours!").c_str());
			}
			size_t found_id = pvec_neighbours_sorted[pid - 1]->p_center_pa->id;
			if (found_id != pid)
			{
				throw std::runtime_error("pid not match when get neighbours through sorted method, this could due to inconsecutive pid");
			}
#endif // SAFE_GET_NEIGHBOURS
			return *pvec_neighbours_sorted[pid - 1];
		}
	}

	VectorDd Configuration_neighbours::to_box_coordination(const Particle& pa) const
	{
		VectorDd pa_vector_cartesian;
		pa_vector_cartesian(0) = pa.rx;
		pa_vector_cartesian(1) = pa.ry;
		pa_vector_cartesian(2) = pa.rz;
		VectorDd pa_vector_box = m_vector_cartesian_to_box * pa_vector_cartesian;
		return pa_vector_box;
	}

	VectorDd Configuration_neighbours::get_PBC_rij(const Particle& pai, const Particle& paj) const
	{
		VectorDd pai_vector_box = to_box_coordination(pai);
		VectorDd paj_vector_box = to_box_coordination(paj);
		VectorDd rij_box = pai_vector_box - paj_vector_box;

		for (size_t i = 0; i < static_cast<size_t>(rij_box.size()); i++)
		{
			rij_box[i] -= round(rij_box[i]);
		}
		;
		VectorDd rij_cartesian = m_vector_box_to_cartesian * rij_box;
		return rij_cartesian;
	}


	void D2::compute_Xij()
	{
		const Particle& (Configuration::Configuration:: * get_pa)(size_t) const; //member function pointer, choose diffenrent method of find partcile in config,. according to get_flag_sorted() 
		if (config_t_->get_flag_sorted() && config_t->get_flag_sorted()) // if config has been sorted, then use get_particle_sorted to find particle, this method is faster when call-times are numerous
		{
			get_pa = &Configuration::Configuration::get_particle_sorted;
		}
		else // if config has not been sorted, then use travesal method to find particle
		{
			get_pa = &Configuration::Configuration::get_particle;
		}

		Xij.setZero();
		const Particle& pa_0_t = (config_t->*get_pa) (p_pa->id);
		auto& pa_0_t_neighbours = config_t->get_neighbours(pa_0_t);
		const Particle& pa_0_t_ = (config_t_->*get_pa)(pa_0_t.id);
		for (size_t n = 0; n < pa_0_t_neighbours.vec_neighbours.size(); n++)
		{
			const Particle& pa_n_t = *pa_0_t_neighbours.vec_neighbours[n].p_neighbour;
			const Particle& pa_n_t_ = (config_t_->*get_pa)(pa_n_t.id);
#ifndef NO_STORE_NEIGHBOUR_RIJ // if rij is pre-stored, use it
			const VectorDd& r0n_t = pa_0_t_neighbours.vec_neighbours[n].rij;
#else // if rij is NOT pre-stored, calculate it 
			const VectorDd r0n_t = config_t->get_PBC_rij(pa_0_t, pa_n_t);
#endif // !NO_STORE_NEIGHBOUR_RIJ
			const VectorDd r0n_t_ = config_t_->get_PBC_rij(pa_0_t_, pa_n_t_);
			Xij += r0n_t * r0n_t_.transpose();
		}
		flag_computed_Xij = true;
		return;
	}

	void D2::compute_Yij()
	{
		const Particle& (Configuration::Configuration:: * get_pa)(size_t) const;
		if (config_t_->get_flag_sorted() && config_t->get_flag_sorted())
		{
			get_pa = &Configuration::Configuration::get_particle_sorted;
		}
		else
		{
			get_pa = &Configuration::Configuration::get_particle;
		}

		Yij.setZero();
		const Particle& pa_0_t = (config_t->*get_pa) (p_pa->id);
		auto& pa_0_t_neighbours = config_t->get_neighbours(pa_0_t);

		const Particle& pa_0_t_ = (config_t_->*get_pa)(pa_0_t.id);
		for (size_t n = 0; n < pa_0_t_neighbours.vec_neighbours.size(); n++)
		{
			const Particle& pa_n_t = *pa_0_t_neighbours.vec_neighbours[n].p_neighbour;
			const Particle& pa_n_t_ = (config_t_->*get_pa)(pa_n_t.id);
			const VectorDd rn0_t_ = config_t_->get_PBC_rij(pa_n_t_, pa_0_t_);
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
		const Particle& (Configuration::Configuration:: * get_pa)(size_t) const;
		if (config_t_->get_flag_sorted() && config_t->get_flag_sorted())
		{
			get_pa = &Configuration::Configuration::get_particle_sorted;
		}
		else
		{
			get_pa = &Configuration::Configuration::get_particle;
		}

		compute_epsilonij();
		const Particle& pa_0_t = (config_t->*get_pa) (p_pa->id);
		const Particle& pa_0_t_ = (config_t_->*get_pa)(pa_0_t.id);
		const auto& pvec_ngbrs_0_t = config_t->get_neighbours(pa_0_t).vec_neighbours;
		d2 = 0;
		for (size_t n = 0; n < pvec_ngbrs_0_t.size(); n++)
		{
			const Particle& pa_n_t = *pvec_ngbrs_0_t[n].p_neighbour;
			const Particle& pa_n_t_ = (config_t_->*get_pa)(pa_n_t.id);
#ifndef NO_STORE_NEIGHBOUR_RIJ // if rij is pre-stored, use it
			const VectorDd& r0n_t = pvec_ngbrs_0_t[n].rij;
#else // if rij is NOT pre-stored, calculate it
			const VectorDd r0n_t = config_t->get_PBC_rij(pa_0_t, pa_n_t);
#endif // !NO_STORE_NEIGHBOUR_RIJ
			const VectorDd r0n_t_ = config_t_->get_PBC_rij(pa_0_t_, pa_n_t_);
			VectorDd Dn_vector = r0n_t - (DELTAij + epsilonij) * r0n_t_;
			double Dn = Dn_vector.dot(Dn_vector);
			d2 += Dn;
		}
		d2 /= pvec_ngbrs_0_t.size();
		flag_computed_d2 = true;
	}
}
