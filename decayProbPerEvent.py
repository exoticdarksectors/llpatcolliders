import numpy as np
import pandas as pd
import trimesh
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from tqdm import tqdm
from scipy.integrate import quad
from plottingFunctions import plot_tube_shape,plot_tube_interactive


# Speed of light in m/s
SPEED_OF_LIGHT = 299792458.0  # m/s

def eta_phi_to_direction(eta, phi):
    """
    Convert pseudorapidity (eta) and azimuthal angle (phi) to 3D direction vector
    
    Args:
        eta: Pseudorapidity
        phi: Azimuthal angle in radians
    
    Returns:
        Normalized 3D direction vector [dx, dy, dz]
    """
    # Convert eta to theta (polar angle)
    theta = 2 * np.arctan(np.exp(-eta))
    
    # Convert to Cartesian direction
    dx = np.sin(theta) * np.cos(phi)
    dy = np.sin(theta) * np.sin(phi)
    dz = np.cos(theta)
    
    return np.array([dx, dy, dz])

def calculate_decay_length(momentum, mass, lifetime):
    """
    Calculate decay length for a particle
    
    Args:
        momentum: Particle momentum in GeV/c
        mass: Particle mass in GeV/c²
        lifetime: Particle lifetime in seconds
    
    Returns:
        Decay length in meters
    """
    # Calculate energy
    energy = np.sqrt(momentum**2 + mass**2)
    
    # Calculate velocity (as fraction of c)
    beta = momentum / energy
    
    # Calculate gamma factor
    gamma = energy / mass
    
    # Decay length = gamma * beta * c * tau
    decay_length = gamma * beta * SPEED_OF_LIGHT * lifetime
    
    return decay_length

def process_particle_csv(csv_file, mesh, origin, lifetime_seconds):
    """
    Process CSV file containing particle data and calculate decay probabilities
    
    Args:
        csv_file: Path to CSV file with columns: event, eta, phi, momentum, mass
        mesh: Trimesh object representing the tube
        origin: Origin point [x, y, z]
        lifetime_seconds: Particle lifetime in seconds
    
    Returns:
        DataFrame with original data plus decay probabilities and event-level statistics
    """
    # Read CSV file
    df = pd.read_csv(csv_file)
    
    # Ensure required columns exist
    required_columns = ['event', 'eta', 'phi', 'momentum', 'mass']
    if not all(col in df.columns for col in required_columns):
        raise ValueError(f"CSV must contain columns: {required_columns}")
    
    # Initialize particle-level results columns
    df['decay_length'] = np.zeros(len(df))
    df['hits_tube'] = False
    df['entry_distance'] = np.nan
    df['exit_distance'] = np.nan
    df['path_length_in_tube'] = np.nan
    df['decay_probability'] = 0.0
    
    origin = np.array(origin)
    
    print(f"Processing {len(df)} particles from {df['event'].nunique()} events...")
    
    # Process each particle
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing particles"):
        # Get particle properties
        eta = row['eta']
        phi = row['phi']
        momentum = row['momentum']  # GeV/c
        mass = row['mass']  # GeV/c²
        
        # Calculate decay length
        decay_length = calculate_decay_length(momentum, mass, lifetime_seconds)
        df.at[idx, 'decay_length'] = decay_length
        
        # Convert eta, phi to direction vector
        direction = eta_phi_to_direction(eta, phi)
        
        # Find intersections with tube
        locations, _, _ = mesh.ray.intersects_location(
            ray_origins=[origin],
            ray_directions=[direction]
        )
        
        if len(locations) >= 2:
            # Particle hits the tube
            df.at[idx, 'hits_tube'] = True
            
            # Calculate distances
            distances = [np.linalg.norm(loc - origin) for loc in locations]
            sorted_distances = sorted(distances)
            
            entry_distance = sorted_distances[0]
            exit_distance = sorted_distances[1]
            path_length = exit_distance - entry_distance
            
            df.at[idx, 'entry_distance'] = entry_distance
            df.at[idx, 'exit_distance'] = exit_distance
            df.at[idx, 'path_length_in_tube'] = path_length
            
            # Calculate decay probability
            # P(decay in tube) = P(survive to entry) * P(decay in tube)
            p_survive_to_entry = np.exp(-entry_distance / decay_length)
            p_decay_in_tube = p_survive_to_entry * (1 - np.exp(-path_length / decay_length))
            
            df.at[idx, 'decay_probability'] = p_decay_in_tube
    
    # Calculate event-level statistics
    event_stats = calculate_event_statistics(df)
    
    # Merge event statistics back to particle dataframe
    df = df.merge(event_stats, on='event', how='left')
    
    return df

def calculate_event_statistics(df):
    """
    Calculate event-level statistics for decay probabilities
    
    Args:
        df: DataFrame with particle-level data
    
    Returns:
        DataFrame with event-level statistics
    """
    event_groups = df.groupby('event')
    
    event_stats = []
    
    for event_idx, group in event_groups:
        n_particles = len(group)
        
        # Get decay probabilities for particles in this event
        p1 = group.iloc[0]['decay_probability'] if n_particles >= 1 else 0
        p2 = group.iloc[1]['decay_probability'] if n_particles >= 2 else 0
        
        # Probability that at least one particle decays
        # P(at least one) = 1 - P(none decay) = 1 - (1-p1)(1-p2)
        p_at_least_one = 1 - (1 - p1) * (1 - p2)
        
        # Probability that both decay
        p_both = p1 * p2
        
        # Probability that exactly one decays
        p_exactly_one = p_at_least_one - p_both
        
        event_stats.append({
            'event': event_idx,
            'n_particles_in_event': n_particles,
            'n_particles_hitting_tube': group['hits_tube'].sum(),
            'prob_at_least_one_decays': p_at_least_one,
            'prob_both_decay': p_both,
            'prob_exactly_one_decays': p_exactly_one,
            'particle1_decay_prob': p1,
            'particle2_decay_prob': p2
        })
    
    return pd.DataFrame(event_stats)

def analyze_decay_vs_lifetime(csv_file, mesh, origin, lifetime_range):
    """
    Analyze how decay probability varies with lifetime for particles in CSV
    
    Args:
        csv_file: Path to CSV file
        mesh: Trimesh object
        origin: Origin point
        lifetime_range: Array of lifetimes to test (in seconds)
    
    Returns:
        Dictionary with analysis results
    """
    # Read CSV once to get basic info
    df_base = pd.read_csv(csv_file)
    n_events = df_base['event'].nunique()
    
    results = {
        'lifetimes': lifetime_range,
        'mean_single_particle_decay_prob': [],
        'mean_at_least_one_decay_prob': [],
        'mean_both_decay_prob': [],
        'frac_events_with_decay': [],
        'exclusion': [],
        'total_events': n_events
    }
    
    # For each lifetime, calculate decay probabilities
    for lifetime in tqdm(lifetime_range, desc="Scanning lifetimes"):
        df = process_particle_csv(csv_file, mesh, origin, lifetime)
        
        # Get event statistics
        event_stats = df.groupby('event').first()
        
        # Calculate statistics
        mean_single = df[df['hits_tube']]['decay_probability'].mean() if any(df['hits_tube']) else 0
        mean_at_least_one = event_stats['prob_at_least_one_decays'].mean()
        mean_both = event_stats['prob_both_decay'].mean()
        frac_with_decay = (event_stats['prob_at_least_one_decays'] > 0.01).mean()  # Events with >1% decay prob
        
        results['mean_single_particle_decay_prob'].append(mean_single)
        results['mean_at_least_one_decay_prob'].append(mean_at_least_one)
        results['exclusion'].append(3/(mean_at_least_one*3000*52E3))
        results['mean_both_decay_prob'].append(mean_both)
        results['frac_events_with_decay'].append(frac_with_decay)
    
    return results

def create_sample_csv(filename, n_events=500):
    """
    Create a sample CSV file with particle data for testing
    Two particles per event
    
    Args:
        filename: Output filename
        n_events: Number of events to generate (2 particles per event)
    """
    np.random.seed(42)  # For reproducibility
    
    data = {
        'event': [],
        'eta': [],
        'phi': [],
        'momentum': [],
        'mass': []
    }
    
    for event_idx in range(n_events):
        # Generate two particles per event
        for particle in range(2):
            data['event'].append(event_idx)
            
            # Generate correlated eta values (particles from same decay often go in similar directions)
            if particle == 0:
                base_eta = np.random.uniform(-2.5, 2.5)
                base_phi = np.random.uniform(0, 2*np.pi)
            else:
                # Second particle somewhat correlated with first
                data['eta'].append(base_eta + np.random.normal(0, 0.5))
                data['phi'].append((base_phi + np.pi + np.random.normal(0, 0.3)) % (2*np.pi))
                continue
            
            data['eta'].append(base_eta)
            data['phi'].append(base_phi)
        
        # Momentum and mass (same for both particles in event for simplicity)
        momentum = np.random.lognormal(np.log(10), 0.5)
        mass = np.random.choice([0.140, 0.494, 0.938, 1.116], p=[0.6, 0.2, 0.15, 0.05])
        
        data['momentum'].extend([momentum, momentum])
        data['mass'].extend([mass, mass])
    
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)
    print(f"Created sample CSV with {n_events} events ({len(df)} particles): {filename}")
    return df

# Create the tube mesh (from your provided code)
def create_tube_mesh(path_points, radius=1.0, n_segments=16):
    """Create a tube mesh along a path with circular cross-section"""
    vertices = []
    faces = []
    
    for i in range(len(path_points)):
        if i == 0:
            tangent = path_points[1] - path_points[0]
        elif i == len(path_points) - 1:
            tangent = path_points[i] - path_points[i-1]
        else:
            tangent = path_points[i+1] - path_points[i-1]
        
        tangent = tangent / np.linalg.norm(tangent)
        
        if abs(tangent[2]) < 0.9:
            up = np.array([0, 0, 1])
        else:
            up = np.array([1, 0, 0])
        
        right = np.cross(tangent, up)
        right = right / np.linalg.norm(right)
        up = np.cross(right, tangent)
        up = up / np.linalg.norm(up)
        
        for j in range(n_segments):
            angle = 2 * np.pi * j / n_segments
            offset = radius * (np.cos(angle) * right + np.sin(angle) * up)
            vertex = path_points[i] + offset
            vertices.append(vertex)
        
        if i > 0:
            for j in range(n_segments):
                v1 = (i-1) * n_segments + j
                v2 = (i-1) * n_segments + (j + 1) % n_segments
                v3 = i * n_segments + (j + 1) % n_segments
                v4 = i * n_segments + j
                
                faces.append([v1, v4, v3])
                faces.append([v1, v3, v2])
    
    # Cap the ends
    center_start = len(vertices)
    vertices.append(path_points[0])
    for j in range(n_segments):
        v1 = j
        v2 = (j + 1) % n_segments
        faces.append([center_start, v1, v2])
    
    center_end = len(vertices)
    vertices.append(path_points[-1])
    last_ring_start = (len(path_points) - 1) * n_segments
    for j in range(n_segments):
        v1 = last_ring_start + j
        v2 = last_ring_start + (j + 1) % n_segments
        faces.append([center_end, v2, v1])
    
    return np.array(vertices), np.array(faces)

# Setup tube geometry (from your code)
correctedVert = [
(-86.57954338701529, 0.1882163986665546  ),
(-1731.590867740335, 3.764327973349282   ),
(-3549.761278867689, 7.716872345365118   ),
(-5887.408950317142, 12.798715109387558  ),
(-8053.403266181902, -504.23173203003535 ),
(-10046.991360867298, -1282.5065405198511),
(-11783.350377373874, -2930.9057600491833),
(-12913.652590171332, -4580.622494369192 ),
(-13095.344153684957, -7536.749251839814 ),
(-13099.610392054752, -9015.000846973791 ),
(-13278.792403586143, -11101.567842600896),
(-13372.39869252341, -13536.146959364076 ),
(-13292.093029091975, -15710.234580371536),
(-12779.140603923677, -17972.21925955668   ), 
(-11659.12755425337, -19887.69754879509    ),
(-10105.714877251532, -21630.204967658145  ),
(-7512.845769209047, -23201.0590309365     ),
(-5262.530506741277, -23466.820585854904   ),
(-2751.72374851779, -23472.278861416264    ),
(-241.41890069074725, -23651.64908934632   ),
(1749.6596420124115, -23742.93404270002    ),
(3827.568683300815, -23747.45123626804     ),
(6078.6368113632525, -23752.344862633392   ),
(8502.613071001502, -23844.570897980426    ),
(11446.568501358292, -23764.01427935077    ),
(13438.399909656131, -23594.431304151418   ),
(15777.051401898476, -23251.689242178036   ),
(18289.614846509525, -22648.455684448927   ),
(20889.761655300477, -21697.58643838109    ),
(23143.841245741598, -20659.00835053422    ),
(25486.006110759066, -19098.88262197991    ),
(27742.09334278597, -17364.656724658227    ),
(28871.391734790544, -16062.763895075637   ),
(30781.662703665817, -14153.873179790575   ),
(32518.021720172394, -12505.473960261239   ),
(34513.49197884447, -11075.029330388788    ),
(36636.57295581305, -10427.47081077351     ),
(38759.40297758341, -9866.868267342572     ),
(41357.416667189485, -9655.12481884172     ),
(43694.93886103982, -9703.684649697909     ),
(46379.03018363646, -9666.041369964427     ),
(49409.43967978114, -9629.150955825604     ),
(51660.88424064092, -9503.610617914434     ),
(54258.0195870532, -9596.213086058811      ),
(57028.564975437745, -9602.236010816167    ),
(59539.87364405768, -9433.782334008818     ),
(62050.42944708294, -9526.196585754526     )]

correctedVertWithShift = []
for x,y in correctedVert:
    correctedVertWithShift.append(((x- 11908.8279764855)/1000,(y+13591.106147774964)/1000))

# Create 3D path at z=22
path_vertices = correctedVertWithShift
Z_POSITION = 22
path_3d = np.array([[x, y, Z_POSITION] for x, y in path_vertices])

# Create the tube mesh
tube_radius = 1.4*1.1
vertices, faces = create_tube_mesh(path_3d, radius=tube_radius, n_segments=32)
mesh = trimesh.Trimesh(vertices=vertices, faces=faces)

if mesh.volume < 0:
    mesh.invert()

print("Tube mesh created at z=22")
print(f"Mesh bounds: X:[{mesh.bounds[0][0]:.1f}, {mesh.bounds[1][0]:.1f}], "
      f"Y:[{mesh.bounds[0][1]:.1f}, {mesh.bounds[1][1]:.1f}], "
      f"Z:[{mesh.bounds[0][2]:.1f}, {mesh.bounds[1][2]:.1f}]")

# Main analysis
if __name__ == "__main__":
    # Create sample CSV if needed
    sample_csv = "LLP.csv"
    # create_sample_csv(sample_csv, n_events=500)
    # Plot the tube
    
    # Set origin
    origin = [0, 0, 0]
    plot_tube_shape(mesh, path_3d, tube_radius, origin,Z_POSITION)
    # exit()
    
    # Single lifetime analysis
    print("\n" + "="*50)
    print("SINGLE LIFETIME ANALYSIS")
    print("="*50)
    
    lifetime = 1e-8  # 10 nanoseconds
    df_results = process_particle_csv(sample_csv, mesh, origin, lifetime)
    
    # Print summary statistics
    event_stats = df_results.groupby('event').first()
    
    print(f"\nResults for lifetime = {lifetime*1e9:.1f} nanoseconds:")
    print(f"Total events: {len(event_stats)}")
    print(f"Total particles: {len(df_results)}")
    
    particles_hitting = df_results[df_results['hits_tube'] == True]
    print(f"\nParticle-level statistics:")
    print(f"Particles hitting tube: {len(particles_hitting)} ({len(particles_hitting)/len(df_results)*100:.1f}%)")
    if len(particles_hitting) > 0:
        print(f"Mean decay probability (for hits): {particles_hitting['decay_probability'].mean():.4f}")
    
    print(f"\nEvent-level statistics:")
    print(f"Events with at least one particle hitting tube: {(event_stats['n_particles_hitting_tube'] > 0).sum()}")
    print(f"Events with both particles hitting tube: {(event_stats['n_particles_hitting_tube'] == 2).sum()}")
    print(f"Mean P(at least one decays): {event_stats['prob_at_least_one_decays'].mean():.6f}")
    print(f"Mean P(both decay): {event_stats['prob_both_decay'].mean():.6f}")
    print(f"Mean P(exactly one decays): {event_stats['prob_exactly_one_decays'].mean():.6f}")
    
    # Save results
    df_results.to_csv("particle_decay_results_with_events.csv", index=False)
    event_stats.to_csv("event_decay_statistics.csv")
    print("\nDetailed results saved to CSV files")
    
    # Lifetime scan
    print("\n" + "="*50)
    print("LIFETIME SCAN ANALYSIS")
    print("="*50)
    
    # Range of lifetimes (in seconds)
    lifetimes = np.logspace(-9.5, -4.5, 20)  # 0.1 ns to 1 microsecond
    
    scan_results = analyze_decay_vs_lifetime(sample_csv, mesh, origin, lifetimes)
    
    # Visualization
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Single particle decay probability
    ax1 = axes[0, 0]
    ax1.semilogx(lifetimes * 1e9, scan_results['mean_single_particle_decay_prob'], 
                 'b-', linewidth=2, label='Single particle (hits only)')
    ax1.set_xlabel('Lifetime (nanoseconds)')
    ax1.set_ylabel('Mean Decay Probability')
    ax1.set_title('Single Particle Decay Probability vs Lifetime')
    ax1.grid(True, which="both", ls="-", alpha=0.2)
    ax1.legend()
    
    # Plot 2: Event-level probabilities
    ax2 = axes[0, 1]
    ax2.loglog(lifetimes * 1e9, scan_results['mean_at_least_one_decay_prob'], 
               'r-', linewidth=2, label='At least one decays')
    ax2.loglog(lifetimes * 1e9, scan_results['mean_both_decay_prob'], 
               'g--', linewidth=2, label='Both decay')
    ax2.set_xlabel('Lifetime (nanoseconds)')
    ax2.set_ylabel('Event Decay Probability')
    ax2.set_title('Event-Level Decay Probabilities vs Lifetime')
    ax2.grid(True, which="both", ls="-", alpha=0.2)
    ax2.legend()
    
    # Plot 3: Comparison of single vs at-least-one
    ax3 = axes[1, 0]
    ax3.loglog(lifetimes * 1e9, scan_results['mean_at_least_one_decay_prob'], 
               'r-', linewidth=2, label='P(≥1 decay per event)')
    ax3.loglog(lifetimes * 1e9, scan_results['mean_single_particle_decay_prob'], 
               'b--', linewidth=2, label='P(single particle decay)')
    ax3.set_xlabel('Lifetime (nanoseconds)')
    ax3.set_ylabel('Decay Probability')
    ax3.set_title('Event vs Single Particle Decay Probabilities')
    ax3.grid(True, which="both", ls="-", alpha=0.2)
    ax3.legend()
    
    # Plot 4: Exclusion
    ax4 = axes[1, 1]
    ax4.loglog(lifetimes * 3E8, scan_results['exclusion'], 
                 color='blue', linewidth=2,label="mQ")
    ax4.set_xlabel('Lifetime (m)')
    ax4.set_ylabel('BR')
    # ax4.set_title('Fraction of Events with >1% Decay Probability')
    ax4.grid(True, which="both", ls="-", alpha=0.2)
    
    plt.tight_layout()
    import numpy as np
    externalLines = {}
    externalLines["MATHUSLA"] = np.loadtxt("external/MATHUSLA.csv",delimiter=",")
    externalLines["CODEX"] = np.loadtxt("external/CODEX.csv",delimiter=",")
    externalLines["ANUBIS"] = np.loadtxt("external/ANUBIS.csv",delimiter=",")
    externalLines["ANUBISOpt"] = np.loadtxt("external/ANUBISOpt.csv",delimiter=",")
    externalLines["ANUBISCons"] = np.loadtxt("external/ANUBISUpdateCons.csv",delimiter=",")
    ax4.loglog(externalLines["MATHUSLA"][:,0],externalLines["MATHUSLA"][:,1],
                 color="green", linewidth=2,label="MATHUSLA")
    ax4.loglog(externalLines["CODEX"][:,0],externalLines["CODEX"][:,1],
                 color="cyan", linewidth=2,label="CODEX b")
    ax4.loglog(externalLines["ANUBIS"][:,0],externalLines["ANUBIS"][:,1],
                 color="purple", linewidth=2,label="ANUBIS")
    ax4.loglog(externalLines["ANUBISOpt"][:,0],externalLines["ANUBISOpt"][:,1],
                 color="purple", linewidth=2,linestyle="--",label="ANUBISOpt")
    ax4.loglog(externalLines["ANUBISCons"][:,0],externalLines["ANUBISCons"][:,1],
                 color="magenta", linewidth=2,linestyle="--",label="ANUBISCons")
    plt.legend()

    plt.savefig('event_exclusion_vs_lifetime.png', dpi=150)
    plt.show()
    
    # Create event visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    exit()
    
    # Select a specific lifetime for visualization
    viz_lifetime = 1e-8  # 10 ns
    df_viz = process_particle_csv(sample_csv, mesh, origin, viz_lifetime)
    
    # Plot 1: Event decay probabilities histogram
    event_probs = df_viz.groupby('event')['prob_at_least_one_decays'].first()
    ax1.hist(event_probs, bins=50, edgecolor='black', alpha=0.7)
    ax1.set_xlabel('P(at least one decay)')
    ax1.set_ylabel('Number of Events')
    ax1.set_title(f'Distribution of Event Decay Probabilities (τ = {viz_lifetime*1e9:.1f} ns)')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Correlation between two particles
    events_both_hit = df_viz[df_viz['n_particles_hitting_tube'] == 2].groupby('event').first()
    if len(events_both_hit) > 0:
        ax2.scatter(events_both_hit['particle1_decay_prob'], 
                   events_both_hit['particle2_decay_prob'], 
                   alpha=0.5, s=20)
        ax2.set_xlabel('Particle 1 Decay Probability')
        ax2.set_ylabel('Particle 2 Decay Probability')
        ax2.set_title('Correlation of Decay Probabilities within Events')
        ax2.grid(True, alpha=0.3)
        
        # Add diagonal line
        max_prob = max(events_both_hit['particle1_decay_prob'].max(), 
                      events_both_hit['particle2_decay_prob'].max())
        ax2.plot([0, max_prob], [0, max_prob], 'r--', alpha=0.5, label='Equal probability')
        ax2.legend()
    
    plt.tight_layout()
    plt.savefig('event_correlation_analysis.png', dpi=150)
    plt.show()
    
    # Summary statistics as a function of lifetime
    print("\n" + "="*40)
    print("SUMMARY: Event Decay Probabilities")
    print("="*40)
    print(f"{'Lifetime (ns)':>15} | {'Mean P(≥1)':>12} | {'Mean P(both)':>12} | {'Events >1%':>12}")
    print("-"*60)
    
    for i in range(0, len(lifetimes), 2):  # Print every other value for brevity
        lt_ns = lifetimes[i] * 1e9
        p_one = scan_results['mean_at_least_one_decay_prob'][i]
        p_both = scan_results['mean_both_decay_prob'][i]
        frac = scan_results['frac_events_with_decay'][i]
        print(f"{lt_ns:>15.2f} | {p_one:>12.6f} | {p_both:>12.6f} | {frac:>12.3f}")
