#include "./mechanics.h"

void epithelial_special_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	// Use this for homogenous attachment
	// BM adhesion 
		// is it time to detach (attachment lifetime)
		// am I unattached by capable? 
			// search through neighbors, find closest BM type agent 
			// form adhesion 
		// elastic adhesion 
	
	// plasto-elastic. 
		// elastic: movement towards rest position 
	
	static int nRP = 0; // "rest_position"
	
	std::vector<double> displacement = pCell->custom_data.vector_variables[nRP].value ; 
	// std::cout << "in special mechanics method for cell id " << pCell->ID << std::endl;
	// calling heterotypic cell update
	// heterotypic_update_cell_velocity(pCell, phenotype, dt);

	// trying homogenous update
	// basically in this, get all the ones which are of same type (epithelial),
	// have adhesion force and repulsive force even between them,  

	custom_cell_update_mechanics( pCell , phenotype , dt );
	
	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}
	
	pCell->state.simple_pressure = 0.0; 
	
	//First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		// pCell->add_potentials(*neighbor);
		// for each neighbor in current voxel, add it's contribution to displacement and velocity
		// add_heterotypic_potentials( pCell, *neighbor ); 
		add_spring_potentials(pCell, *neighbor );
	}
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
		// search for neighboring voxels, if the cell is on a voxel phase and in contact with neighboring voxel subcell
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
			// pCell->add_potentials(*neighbor);
			// for each cell in neigboring voxel, add it's contribution to displacement and velocity
			// add_heterotypic_potentials( pCell, *neighbor ); 
			add_spring_potentials(pCell, *neighbor );
		}
	}

	pCell->update_motility_vector(dt); 
	pCell->velocity += phenotype.motility.motility_vector; 
	
	return; 
}

void add_spring_potentials(Cell* my_cell, Cell* other_agent)
{
	// basically check if the adjacent cell is of same type and connected to link
	// if so, then put a repulsive force between these 2, for now use same repulsion co efficient 
	
	static int nCellID = my_cell->custom_data.find_variable_index( "cell_ID" ); 
	
	// if( this->ID == other_agent->ID )
	// if same cell, return
	if( my_cell == other_agent )
	{ return; }
	
	static int nRAOC = my_cell->custom_data.find_variable_index( "relative_adhesion_other_cells" ); 
	static int nRAOCT = my_cell->custom_data.find_variable_index( "relative_adhesion_other_cell_types" ); 

	static int subcell_repulsion_factor_index = my_cell->custom_data.find_variable_index("subcell_repulsion_factor");
	static int subcell_adhesion_factor_index = my_cell->custom_data.find_variable_index("subcell_adhesion_factor");

	double subcell_repulsion_factor = my_cell->custom_data[subcell_repulsion_factor_index];
	double subcell_adhesion_factor = my_cell->custom_data[subcell_adhesion_factor_index];
	
	double rel_heterotypic_adhesion = my_cell->custom_data[nRAOCT];  
	double rel_other_cells_adhesion = my_cell->custom_data[nRAOC];  
	
	// 12 uniform neighbors at a close packing distance, after dividing out all constants
	static double simple_pressure_scale = 0.027288820670331; // 12 * (1 - sqrt(pi/(2*sqrt(3))))^2 
	// 9.820170012151277; // 12 * ( 1 - sqrt(2*pi/sqrt(3)))^2

	double distance = 0; 
	for( int i = 0 ; i < 3 ; i++ ) 
	{ 
		// find distance between cells 
		my_cell->displacement[i] = my_cell->position[i] - (*other_agent).position[i]; 
		distance += (my_cell->displacement[i]) * (my_cell->displacement[i]); 
	}
	// Make sure that the distance is not zero
	
	distance = std::max(sqrt(distance), 0.00001); 
	
	//Repulsive
	// Repulsion attributes
	double R = my_cell->phenotype.geometry.radius + (*other_agent).phenotype.geometry.radius; 
	
	double RN = my_cell->phenotype.geometry.nuclear_radius + (*other_agent).phenotype.geometry.nuclear_radius;	
	double temp_r, c;
	if( distance > R ) 
	// R is basically distance between centers when they are attached, if distance is more, 
	// then there's still space between them
	{
		temp_r=0;
	}
	else
	{
		temp_r = -distance; // -d
		temp_r /= R; // -d/R
		temp_r += 1.0; // 1-d/R
		temp_r *= temp_r; // (1-d/R)^2 
		temp_r *= subcell_repulsion_factor;
		// temp_r *= 1.5
		
		// add the relative pressure contribution 
		my_cell->state.simple_pressure += ( temp_r / simple_pressure_scale ); // New July 2017 
	}
	
	// August 2017 - back to the original if both have same coefficient 
	double effective_repulsion = sqrt( my_cell->phenotype.mechanics.cell_cell_repulsion_strength * other_agent->phenotype.mechanics.cell_cell_repulsion_strength );
	temp_r *= effective_repulsion; 
	
	//////////////////////////////////////////////////////////////////
	
	// Adhesive
	double max_interactive_distance = my_cell->phenotype.mechanics.relative_maximum_adhesion_distance * my_cell->phenotype.geometry.radius + 
		(*other_agent).phenotype.mechanics.relative_maximum_adhesion_distance * (*other_agent).phenotype.geometry.radius;
		
	if(distance < max_interactive_distance ) 
	{	
		// double temp_a = 1 - distance/max_interactive_distance; 
		double temp_a = -distance; // -d
		temp_a /= max_interactive_distance; // -d/S
		temp_a += 1.0; // 1 - d/S 
		temp_a *= temp_a; // (1-d/S)^2 
		temp_a *= subcell_adhesion_factor;
		
		// August 2017 - back to the original if both have same coefficient 
		double effective_adhesion = sqrt( my_cell->phenotype.mechanics.cell_cell_adhesion_strength * other_agent->phenotype.mechanics.cell_cell_adhesion_strength ); 
		
		int my_id = (int) my_cell->custom_data[nCellID] ; 
		int other_id = (int) other_agent->custom_data[nCellID] ; 
	
		if( my_id != other_id )
		{ effective_adhesion *= rel_other_cells_adhesion; }
		
		if( my_cell->type != other_agent->type )
		{ effective_adhesion *= rel_heterotypic_adhesion; }
		temp_a *= effective_adhesion; 
		// temp_a *= 1.25;
		temp_r -= temp_a;
	}
	/////////////////////////////////////////////////////////////////
	if( fabs(temp_r) < 1e-16 )
	{ return; }
	temp_r /= distance;
	// multiply the temp_r value with velocity and add it to displacement, put it back in displacement
	axpy( &(my_cell->velocity) , temp_r , my_cell->displacement ); 

	static int attach_to_BM_i = my_cell->custom_data.find_variable_index( "attach_to_BM" ); 
	// if different type or same cell, return
	// if (my_cell == other_agent) //|| my_cell->type != other_agent -> type) 
	// {

	// 	return;
	// }
	
	// if the other agent is not spring linked, return
	// if (other_agent->custom_data[attach_to_BM_i] == 0.0)
	// {
	// 	return;
	// }

	// If the other agent is closer to membrane than current cell, then return
	int pbmIndex = microenvironment.find_density_index("pbm");
	int vi_myCell = microenvironment.nearest_voxel_index(my_cell->position);
	int vi_otherCell = microenvironment.nearest_voxel_index(other_agent->position);

	double signed_dist_myCell = microenvironment.density_vector(vi_myCell).at(pbmIndex);
	double signed_dist_otherCell = microenvironment.density_vector(vi_otherCell).at(pbmIndex);

	if ((signed_dist_myCell > 0 ) || (signed_dist_otherCell > 0))
	{
		return;
	}

	// if (std::abs(signed_dist_otherCell) < std::abs(signed_dist_myCell))
	// {
	// 	return;
	// }

	// std::cout << "reached the potential part for cell " << my_cell->ID ;
	// std::cout << " curr sd = " << signed_dist_myCell << " " << "other sd = " << signed_dist_otherCell << std::endl;
	// double distance = 0; 
	// for( int i = 0 ; i < 3 ; i++ ) 
	// { 
	// 	// find distance between cells 
	// 	my_cell->displacement[i] = my_cell->position[i] - (*other_agent).position[i]; 
	// 	distance += (my_cell->displacement[i]) * (my_cell->displacement[i]); 
	// }

	// distance = std::max(sqrt(distance), 0.00001); 
	// double R = my_cell->phenotype.geometry.radius + (*other_agent).phenotype.geometry.radius; 

	// double temp_r;
	// if( distance < R ) 
	// {

		// temp_r = -distance; // -d
		// temp_r /= R; // -d/R
		// temp_r += 1.0; // 1-d/R
		// temp_r *= temp_r; // (1-d/R)^2 
		// // temp_r /= 2;

		// // push the current cell away
		// // disabling force due to linkage
		
		// axpy(&(my_cell->velocity), temp_r, my_cell->displacement);
	// }
	return;
}

void plasto_elastic_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	// BM adhesion 
		// is it time to detach (attachment lifetime)
		// am I unattached by capable? 
			// search through neighbors, find closest BM type agent 
			// form adhesion 
		// elastic adhesion 
	
	// plasto-elastic. 
		// elastic: movement towards rest position 
	
	static int nRP = 0; // "rest_position"
	
	// displacement 
	std::vector<double> displacement = pCell->custom_data.vector_variables[nRP].value - pCell->position; 
	
	static int nEConst = pCell->custom_data.find_variable_index( "cell_elasticity" );
	static int nPConst = pCell->custom_data.find_variable_index( "cell_plasticity" );

	// first, update the agent's velocity based upon the elastic model
	axpy( &( pCell->velocity ) , pCell->custom_data[nEConst] , displacement );

	// now, plastic mechanical relaxation

	double plastic_temp_constant = -dt * pCell->custom_data[nPConst];
	axpy( &(pCell->custom_data.vector_variables[nRP].value) , plastic_temp_constant , displacement );
	
	return; 
}

// specialized potential function 

void add_heterotypic_potentials(Cell* my_cell , Cell* other_agent)
{
	static int nCellID = my_cell->custom_data.find_variable_index( "cell_ID" ); 
	
	// if( this->ID == other_agent->ID )
	if( my_cell == other_agent )
	{ return; }
	
	static int nRAOC = my_cell->custom_data.find_variable_index( "relative_adhesion_other_cells" ); 
	static int nRAOCT = my_cell->custom_data.find_variable_index( "relative_adhesion_other_cell_types" ); 
	
	double rel_heterotypic_adhesion = my_cell->custom_data[nRAOCT];  
	double rel_other_cells_adhesion = my_cell->custom_data[nRAOC];  
	
	// 12 uniform neighbors at a close packing distance, after dividing out all constants
	static double simple_pressure_scale = 0.027288820670331; // 12 * (1 - sqrt(pi/(2*sqrt(3))))^2 
	// 9.820170012151277; // 12 * ( 1 - sqrt(2*pi/sqrt(3)))^2

	double distance = 0; 
	for( int i = 0 ; i < 3 ; i++ ) 
	{ 
		// find distance between cells 
		my_cell->displacement[i] = my_cell->position[i] - (*other_agent).position[i]; 
		distance += (my_cell->displacement[i]) * (my_cell->displacement[i]); 
	}
	// Make sure that the distance is not zero
	
	distance = std::max(sqrt(distance), 0.00001); 
	
	//Repulsive
	// Repulsion attributes
	double R = my_cell->phenotype.geometry.radius + (*other_agent).phenotype.geometry.radius; 
	
	double RN = my_cell->phenotype.geometry.nuclear_radius + (*other_agent).phenotype.geometry.nuclear_radius;	
	double temp_r, c;
	if( distance > R ) 
	// R is basically distance between centers when they are attached, if distance is more, 
	// then there's still space between them
	{
		temp_r=0;
	}
	else
	{
		temp_r = -distance; // -d
		temp_r /= R; // -d/R
		temp_r += 1.0; // 1-d/R
		temp_r *= temp_r; // (1-d/R)^2 
		
		// add the relative pressure contribution 
		my_cell->state.simple_pressure += ( temp_r / simple_pressure_scale ); // New July 2017 
	}
	
	// August 2017 - back to the original if both have same coefficient 
	double effective_repulsion = sqrt( my_cell->phenotype.mechanics.cell_cell_repulsion_strength * other_agent->phenotype.mechanics.cell_cell_repulsion_strength );
	temp_r *= effective_repulsion; 
	
	//////////////////////////////////////////////////////////////////
	
	// Adhesive
	double max_interactive_distance = my_cell->phenotype.mechanics.relative_maximum_adhesion_distance * my_cell->phenotype.geometry.radius + 
		(*other_agent).phenotype.mechanics.relative_maximum_adhesion_distance * (*other_agent).phenotype.geometry.radius;
		
	if(distance < max_interactive_distance ) 
	{	
		// double temp_a = 1 - distance/max_interactive_distance; 
		double temp_a = -distance; // -d
		temp_a /= max_interactive_distance; // -d/S
		temp_a += 1.0; // 1 - d/S 
		temp_a *= temp_a; // (1-d/S)^2 
		
		// August 2017 - back to the original if both have same coefficient 
		double effective_adhesion = sqrt( my_cell->phenotype.mechanics.cell_cell_adhesion_strength * other_agent->phenotype.mechanics.cell_cell_adhesion_strength ); 
		
		int my_id = (int) my_cell->custom_data[nCellID] ; 
		int other_id = (int) other_agent->custom_data[nCellID] ; 
	
		if( my_id != other_id )
		{ effective_adhesion *= rel_other_cells_adhesion; }
		
		if( my_cell->type != other_agent->type )
		{ effective_adhesion *= rel_heterotypic_adhesion; }
		temp_a *= effective_adhesion; 
		
		temp_r -= temp_a;
	}
	/////////////////////////////////////////////////////////////////
	if( fabs(temp_r) < 1e-16 )
	{ return; }
	temp_r /= distance;
	// multiply the temp_r value with velocity and add it to displacement, put it back in displacement
	axpy( &(my_cell->velocity) , temp_r , my_cell->displacement ); 
	
	return;
}

void heterotypic_update_cell_velocity( Cell* pCell, Phenotype& phenotype, double dt)
{
	if( pCell->functions.add_cell_basement_membrane_interactions )
	{
		pCell->functions.add_cell_basement_membrane_interactions(pCell, phenotype,dt);
	}
	
	pCell->state.simple_pressure = 0.0; 
	
	//First check the neighbors in my current voxel
	std::vector<Cell*>::iterator neighbor;
	std::vector<Cell*>::iterator end = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].end();
	for(neighbor = pCell->get_container()->agent_grid[pCell->get_current_mechanics_voxel_index()].begin(); neighbor != end; ++neighbor)
	{
		// pCell->add_potentials(*neighbor);
		// for each neighbor in current voxel, add it's contribution to displacement and velocity
		add_heterotypic_potentials( pCell, *neighbor ); 
	}
	std::vector<int>::iterator neighbor_voxel_index;
	std::vector<int>::iterator neighbor_voxel_index_end = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].end();

	for( neighbor_voxel_index = 
		pCell->get_container()->underlying_mesh.moore_connected_voxel_indices[pCell->get_current_mechanics_voxel_index()].begin();
		neighbor_voxel_index != neighbor_voxel_index_end; 
		++neighbor_voxel_index )
		// search for neighboring voxels, if the cell is on a voxel phase and in contact with neighboring voxel subcell
	{
		if(!is_neighbor_voxel(pCell, pCell->get_container()->underlying_mesh.voxels[pCell->get_current_mechanics_voxel_index()].center, pCell->get_container()->underlying_mesh.voxels[*neighbor_voxel_index].center, *neighbor_voxel_index))
			continue;
		end = pCell->get_container()->agent_grid[*neighbor_voxel_index].end();
		for(neighbor = pCell->get_container()->agent_grid[*neighbor_voxel_index].begin();neighbor != end; ++neighbor)
		{
			// pCell->add_potentials(*neighbor);
			// for each cell in neigboring voxel, add it's contribution to displacement and velocity
			add_heterotypic_potentials( pCell, *neighbor ); 
		}
	}

	pCell->update_motility_vector(dt); 
	pCell->velocity += phenotype.motility.motility_vector; 
	
	return; 
}



void BM_special_mechanics( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	return; 
}

void custom_cell_update_mechanics( Cell* pCell , Phenotype& phenotype , double dt )
{
	static int nRP = 0; // "rest_position"
	
	// displacement 
	std::vector<double> disp = pCell->custom_data.vector_variables[nRP].value - pCell->position; 
	
	static int nEConst = pCell->custom_data.find_variable_index( "cell_elasticity" );
	static int nPConst = pCell->custom_data.find_variable_index( "cell_plasticity" );

	// first, update the agent's velocity based upon the elastic model
	// axpy( &( pCell->velocity ) , pCell->custom_data[nEConst] , disp );

	// now, plastic mechanical relaxation

	// double plastic_temp_constant = -dt * pCell->custom_data[nPConst];
	// axpy( &(pCell->custom_data.vector_variables[nRP].value) , plastic_temp_constant , disp );

	static int membrane_repulsion_factor_index = pCell->custom_data.find_variable_index("membrane_repulsion_factor");
	static int membrane_adhesion_factor_index = pCell->custom_data.find_variable_index("membrane_adhesion_factor");

	double membrane_repulsion_factor = pCell->custom_data[membrane_repulsion_factor_index];
	double membrane_adhesion_factor = pCell->custom_data[membrane_adhesion_factor_index];

    static int attach_lifetime_i = pCell->custom_data.find_variable_index( "attach_lifetime" ); 
    static int attach_time_i = pCell->custom_data.find_variable_index( "attach_time" ); 
    static int attach_to_BM_i = pCell->custom_data.find_variable_index( "attach_to_BM" ); 

	
    // cap letters (X,N, etc) represent vectors
    // Agent vars:
    //   position = X
    //   radius = r
    //   adhesion radius = r_A
    // if (pCell->ID == 4)
    // {
    //     std::cout << " ID = 4, phenotype.geometry.radius = " << phenotype.geometry.radius << std::endl;
    // }

    double adhesion_radius = phenotype.geometry.radius * phenotype.mechanics.relative_maximum_adhesion_distance;
	adhesion_radius = phenotype.geometry.radius * 2;
    int ncells_attached = 0;
	double temp_r;
	double R = pCell->phenotype.geometry.radius * 2;

	int pbmIndex = microenvironment.find_density_index("pbm");
	int n_x_index = microenvironment.find_density_index("n_x");
	int n_y_index = microenvironment.find_density_index("n_y");

	int vi = microenvironment.nearest_voxel_index(pCell->position);

	std::vector<double> nearest_voxel = microenvironment.nearest_density_vector(pCell->position);

	
	// double nx = microenvironment.density_vector(nearest_voxel[0], nearest_voxel[1]).at(n_x_index);
	double nx = microenvironment.density_vector(vi).at(n_x_index);
	double ny = microenvironment.density_vector(vi).at(n_y_index);
	// double ny = microenvironment.density_vector(nearest_voxel[0], nearest_voxel[1]).at(n_y_index);
	double signed_dist = microenvironment.density_vector(vi).at(pbmIndex);
	// double signed_dist = microenvironment.density_vector(nearest_voxel[0], nearest_voxel[1]).at(pbmIndex);

	// double displacement = 0.0 - pCell->position[1];  // displacement: just (negative) y (height) for test case
	double displacement = signed_dist;
	// displacement = std::sqrt((pCell->position[0]- nx) * (pCell->position[0]- nx) + (pCell->position[1]- ny) * (pCell->position[1]- ny));

	double dx = nx - pCell->position[0];
	double dy = ny - pCell->position[1];
	std::vector<double> normal = {dx, dy, 0};
	double dv = 0.01;
	//===================================
    //   attach
    //===================================

	if (PhysiCell_globals.current_time > 6)
	{
		phenotype.motility.is_motile = false;
	}

	if  (pCell->custom_data[attach_to_BM_i] == 1.0 )
	{
		if (displacement < 0) 
		{	
			temp_r = displacement * 0.001; // d
			temp_r /= adhesion_radius; // d/R
			temp_r += 1.0; // 1-d/R 
			temp_r *= temp_r; // (1-d/R)^2 
			temp_r *= membrane_adhesion_factor;
			// temp_r *= -1;
			axpy(&(pCell->velocity), temp_r , normal);
		}
		else 
		{
			// if crosses the barrier, then zoom it back inside
			// std::cout << "crossed the membrane, pushing back in " << pCell->ID << std::endl; 
			pCell->custom_data[attach_to_BM_i] == 0.0;
			pCell->custom_data[attach_time_i] = 0.0;
			dv *= 1000;
			axpy(&(pCell->velocity), dv , normal);
			return;
		}

	}


	if( pCell->custom_data[attach_to_BM_i] == 0.0 )  // not attached to BM
	{
		
        if (displacement <= 0.0 && displacement > -adhesion_radius )
        {
            std::cout << "t="<<PhysiCell_globals.current_time << "attaching ID=" << pCell->ID << ": displacement= " << displacement <<", adhesion radius= " << adhesion_radius << std::endl;
            // double p_BM = pv - d*nv
            pCell->custom_data[attach_to_BM_i] = 1.0;   // attached to BM now
            pCell->custom_data[attach_time_i] = 0.0;   // reset its time of being attached
			std::cout << "velocity " << pCell->velocity[0] << " " << pCell->velocity[1] <<std::endl;

			temp_r = displacement; // d
			temp_r /= adhesion_radius; // d/R
			temp_r += 1.0; // 1-d/R 
			temp_r *= temp_r; // (1-d/R)^2 
			temp_r *= membrane_adhesion_factor;
			// temp_r *= -1;
			

			axpy(&(pCell->velocity), temp_r , normal);
            // phenotype.motility.is_motile = false; 
        }
		// else if (displacement <= 0.0 && displacement < -adhesion_radius ) 
		// {
			
		// }
		
	}

		// double temp_r;
		// double R = pCell->phenotype.geometry.radius * 2;

		temp_r = 1; // 1
		temp_r /= displacement; // 1/d
		temp_r *= membrane_repulsion_factor;
		// push the current cell away
		// disabling force due to linkage
		axpy(&(pCell->velocity), temp_r , normal);


	if (pCell->custom_data[attach_time_i] > pCell->custom_data[attach_lifetime_i])
	{
		pCell->custom_data[attach_to_BM_i] == 0.0;
		pCell->custom_data[attach_time_i] = 0.0;
		return;
	}

	
    
    pCell->custom_data[attach_time_i] += dt;
    
	return; 
}