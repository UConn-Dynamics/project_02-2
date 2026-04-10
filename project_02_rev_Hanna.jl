### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ f5b9dd00-5682-4241-8b80-62ab0ac3bda8
begin
	### A Pluto.jl notebook ###
	# v0.20.25
	
	using Markdown
	using InteractiveUtils
end

# ╔═╡ 5f84fe8f-ee3c-4963-b458-5d5b11450f5d
using LinearAlgebra, Plots, DSP

# ╔═╡ 6ce190b9-7344-4cbf-8e12-adefc8a3bba2
md"""# Project_02 - Multibody kinematic modeling

![Dual slider kinematics project](https://raw.githubusercontent.com/cooperrc/me5180-project_02/refs/heads/main/dual-slider.svg)

In this project, a rigid bar is connected to two sliding pistons along
the diagonal tracks. As the pistons move along the tracks, the rigid bar rotates at a constant rate, $\dot{\theta}_3 = 2~rad/s$. The figure above has three _relative_ ccoordinate systems that move with the bodies:

1. $x_1-y_1-$ describes piston 1 position and orientation, $\theta_1$
2. $x_2-y_2-$ describes piston 2 position and orientation, $\theta_2$
3. $x_3-y_3-$ describes the rigid bar position and orientation, $\theta_3$

Each of the pistons are on tracks at $\pm 45^o$ and the rotating rigid
bar is 10 cm. The hinges are mounted to the center of the pistons
connecting the ends of the rigid bar. 
 
In this project, you need to 

1. determine constraint equations $C(\mathbf{q},~t)$
2. solve for the velocities, $\dot{q}$ and accelerations, $\ddot{q}$
3. visualize the motion of the system as the rigid bar goes through at least one full rotation
"""

# ╔═╡ 943023c9-f07b-486d-9b0e-3983b0147af6
md"""# 1. Determine the Constraint equations C(q,t)


The system consists of two pistons constrained to slide along fixed diagonal tracks and connected by a rigid bar of fixed length L=0.10m
The pistons are connected to the ends of the rigid bar through ideal pin joints located at the piston centers.

The rigid bar rotates with a prescribed constant angular velocity $\theta_3$.

A global inertial coordinate system $(x, y)$ is defined with its origin at the intersection point of the two tracks. Each piston moves along a straight line defined by its track angles $\alpha_1$ and $\alpha_2$, respectively

**Generalized Coordinates:**

Instead of modeling each body with full rigid‑body coordinates, we use reduced generalized coordinates ->

$q=\begin{bmatrix}
s_1 \\
s_2 \\
\theta_3
\end{bmatrix}$

**Geometry:** 

The geometry of $r_1$ and $r_2$ are expressed below

$\mathbf{r}_{1} = s_1\begin{bmatrix}
cos{\theta}_{1} \\
sin{\theta}_{1}
\end{bmatrix}=s_1\begin{bmatrix}
cos(45) \\
sin(45)
\end{bmatrix}$
$\mathbf{r}_{2} = s_2\begin{bmatrix}
cos{\theta}_{2} \\
sin{\theta}_{2}
\end{bmatrix}=s_2\begin{bmatrix}
cos(135) \\
sin(135)
\end{bmatrix}$


**Setting up the Constraint Equation:**

The global position vectors of the pistons are written directly from the track geometry. The rigid bar enforces a distance constraint between the pistons. 

$(\mathbf{r}_{2}-\mathbf{r}_{1})^2=L^2$
$C(q,t)=(\mathbf{r}_{2}-\mathbf{r}_{1})*(\mathbf{r}_{2}-\mathbf{r}_{1})-L^2=0$

$\mathbf{r}_{2}-\mathbf{r}_{1}=\begin{bmatrix}
s_2cos(135)-s_1cos(45) \\
s_2sin(135)-s_1sin(45)
\end{bmatrix}=\begin{bmatrix}
s_2(-\frac{\sqrt{2}}{2})-s_1(\frac{\sqrt{2}}{2}) \\
s_2(\frac{\sqrt{2}}{2})-s_1(\frac{\sqrt{2}}{2})
\end{bmatrix}$


**Final Constraint Equation:**

Substituting the piston position expression gives

$C(q,t)=s_1^2+s_2^2-L^2=0$


**For Julia/Pluto Purposes:**

$C(q,t)=A^2-b^2=0$
$A=\mathbf{r}_{2}-\mathbf{r}_{1}=\begin{bmatrix}
s_2cos(135)-s_1cos(45) \\
s_2sin(135)-s_1sin(45)
\end{bmatrix}$
$b=L=\begin{bmatrix}
Lcos{\theta}_{3} \\
Lsin{\theta}_{3}
\end{bmatrix}$
"""

# ╔═╡ 758462d7-cf01-4f41-9641-b806a3adfbfe
md"""

**Vector from Piston 1 to Piston 2:**


${\theta}_{3}(t)=0$

$\mathbf{r}_{21} = \begin{bmatrix}
Lcos{\theta}_{3} \\
Lsin{\theta}_{3}
\end{bmatrix}$

$\mathbf{r}_{2} = \mathbf{r}_{1}+\mathbf{r}_{21}$

This vector is constrained to remain aligned with the rigid bar and have magnitude L.

**XY Directions:**

The constraint is enforced independently in the x- and y-directions, producing two scaler equations->

$s_2\begin{bmatrix}
cos(135) \\
sin(135)
\end{bmatrix}=s_1\begin{bmatrix}
cos(45) \\
sin(45)
\end{bmatrix} + \begin{bmatrix}
Lcos{\theta}_{3} \\
Lsin{\theta}_{3}
\end{bmatrix}$


$s_2\frac{-2}{\sqrt{2}}=s_1\frac{2}{\sqrt{2}}+Lcos{\theta}_{3}$
$s_2\frac{2}{\sqrt{2}}=s_1\frac{2}{\sqrt{2}}+Lsin{\theta}_{3}$


**Addition Method then isolate s value:**

Solving the coupled linear equations yields the pistion displacements->

s_1:

$0=\frac{4}{\sqrt{2}}s_1+L(cos{\theta}_{3}+sin{\theta}_{3})$
$0=2\sqrt{2}s_1+L(cos{\theta}_{3}+sin{\theta}_{3})$
$s_1=-\frac{L}{2\sqrt{2}}(cos{\theta}_{3}+sin{\theta}_{3})$


s_2:

$\frac{4}{\sqrt{2}}s_2=L(sin{\theta}_{3}-cos{\theta}_{3})$
$2{\sqrt{2}}s_2=L(sin{\theta}_{3}-cos{\theta}_{3})$
$s_2=\frac{L}{2\sqrt{2}}(sin{\theta}_{3}-cos{\theta}_{3})$

"""

# ╔═╡ c46aa017-b129-4596-aa65-15d4d4aca963
md"""

**Acceleration:**

Differentitating once more ->

$C(\ddot{q},t)=\dot{s_1}^2+s_1\ddot{s_1}+\dot{s_2}^2+s_2\ddot{s_2}=0$
$C(\ddot{q},t)=\dot{s_1}^2+\dot{s_2}^2+s_1\ddot{s_1}+s_2\ddot{s_2}=0$

$C(\ddot{q},t)=(\frac{L}{\sqrt{2}}(sin{\theta}_{3}-cos{\theta}_{3}))^2+(\frac{L}{\sqrt{2}}(cos{\theta}_{3}+sin{\theta}_{3}))^2+(-\frac{L}{2\sqrt{2}}(cos{\theta}_{3}+sin{\theta}_{3}))\ddot{s_1}+(\frac{L}{2\sqrt{2}}(sin{\theta}_{3}-cos{\theta}_{3}))\ddot{s_2}=0$

s_1:

$\ddot{s_1}=\frac{2L}{\sqrt{2}}(cos{\theta}_{3}+sin{\theta}_{3})$

s_2:

$\ddot{s_2}=\frac{2L}{\sqrt{2}}(-sin{\theta}_{3}+cos{\theta}_{3})$

"""

# ╔═╡ 8cba0c20-1450-47ce-b2eb-714f18916f4f
# -------------------------
# Parameters
# -------------------------
begin
	L = 0.10                  # bar length in meters (10 cm)
	ω3 = 2.0                  # rad/s
	θ30 = 0.0               # initial angle
	# Verify these track angles from your original assignment
	α1 = π/4                  # example: 45 deg
	α2 = 3*π/4                  # example: 45 deg (shift in angle accounted for later)
end

# ╔═╡ 1b810674-e02c-452e-b32a-d4137b39ea8e
# Time
tspan = range(0, 4π/ω3, length=300)   # one full rotation

# ╔═╡ 694a9ed8-8914-4571-8298-505dd11f6b78
md"""
${\theta}_{3}(t)={\theta}_{3}+{\omega}_{3}t$
"""

# ╔═╡ 370a9b3c-236d-46c7-8d6f-4ee25947ec61
# -------------------------
# Functions
# -------------------------
θ3(t) = θ30 + ω3*t

# ╔═╡ 1e7c0d64-3b10-4162-86c2-9c6b15ec644d
md"""
Position Inputs (Assumes ${\alpha}_{1}$ & ${\alpha_2}$ are equivalent):

$A=\begin{bmatrix}
~-cos({\alpha}_{1})~~~cos({\alpha}_{2}) \\
~-sin({\alpha}_{1})~~sin({\alpha}_{2})
\end{bmatrix}$


$b=\begin{bmatrix}
Lcos({\theta_3}) \\
Lsin({\theta_3})
\end{bmatrix}$

$s=A,~b$
"""

# ╔═╡ 73a6416f-506f-4eec-a03f-a7a9bec789d0
function solve_positions(t, α1, α2, L)
    θ = θ3(t)

    A = [
        -cos(α1)   cos(α2)
        -sin(α1)   sin(α2)
    ]

    b = [
        L*cos(θ)
        L*sin(θ)
    ]

    s = A \ b
    return s[1], s[2]
end

# ╔═╡ 666b9581-f2f9-4c58-b245-59ad60803772
md"""
Velocity Inputs:

$A=\begin{bmatrix}
~-cos({\alpha}_{1})~~~cos({\alpha}_{2}) \\
~-sin({\alpha}_{1})~~sin({\alpha}_{2})
\end{bmatrix}$


$b=\begin{bmatrix}
-L{\omega}_{3}^2sin({\theta_3}) \\
L{\omega}_{3}^2cos({\theta_3})
\end{bmatrix}$

$\dot{s}=A,~b$
"""

# ╔═╡ 179b675e-cad5-47c2-9df1-0c6db8541b79
function solve_velocities(t, α1, α2, L, ω3)
    θ = θ3(t)

    A = [
        -cos(α1)   cos(α2)
        -sin(α1)   sin(α2)
    ]

    b = [
        -L*ω3*sin(θ)
         L*ω3*cos(θ)
    ]

    sdot = A \ b
    return sdot[1], sdot[2]
end

# ╔═╡ ef347f87-3dbc-4fbd-bef5-8008daccfa7a
md"""
Acceleration Inputs:

$A=\begin{bmatrix}
~-cos({\alpha}_{1})~~~cos({\alpha}_{2}) \\
~-sin({\alpha}_{1})~~sin({\alpha}_{2})
\end{bmatrix}$


$b=\begin{bmatrix}
-L{\omega}_{3}^2cos({\theta_3}) \\
-L{\omega}_{3}^2sin({\theta_3})
\end{bmatrix}$

$\ddot{s}=A,~b$

"""

# ╔═╡ a69048f1-fc46-412a-8455-7e7b8048c7ba
function solve_accelerations(t, α1, α2, L, ω3)
    θ = θ3(t)

    A = [
        -cos(α1)   cos(α2)
        -sin(α1)   sin(α2)
    ]

    b = [
        -L*ω3^2*cos(θ)
        -L*ω3^2*sin(θ)
    ]

    sddot = A \ b
    return sddot[1], sddot[2]
end

# ╔═╡ 0e457baa-f346-427e-a5dc-6f2c8be711eb
function center_position(t, L)
    θ = θ3(t)
    x = (L/2)*cos(θ)
    y = (L/2)*sin(θ)
    return x, y
end

# ╔═╡ 72752148-5b7d-454e-b377-0cfaae37b33f
# -------------------------
# Compute results
# -------------------------
s1_vals = Float64[]

# ╔═╡ e8799a11-75a5-4446-9518-d5f2f909c283
s2_vals = Float64[]

# ╔═╡ f8e5cfa9-3765-43e0-9bae-dc245f7c1eff
s1dot_vals = Float64[]

# ╔═╡ 5d2dcd83-64e9-40b1-8e3e-3d3150d06e76
s2dot_vals = Float64[]

# ╔═╡ dcd74d4d-44a8-4ab7-815f-aa7bed0da8d2
s1ddot_vals = Float64[]

# ╔═╡ a8e64e00-b8d2-4cc6-b2c9-148e84364490
s2ddot_vals = Float64[]

# ╔═╡ e3a09adf-a08a-4e38-8245-87d533038b13
md"""# 3. Visualize the Motion of the System as a Rigid Bar Goes through at least one Full Rotation
"""

# ╔═╡ e2836679-2a9e-48fa-8a1b-dd658d07587a
for t in tspan
    s1, s2 = solve_positions(t, α1, α2, L)
    s1dot, s2dot = solve_velocities(t, α1, α2, L, ω3)
    s1ddot, s2ddot = solve_accelerations(t, α1, α2, L, ω3)

    push!(s1_vals, s1)
    push!(s2_vals, s2)
    push!(s1dot_vals, s1dot)
    push!(s2dot_vals, s2dot)
    push!(s1ddot_vals, s1ddot)
    push!(s2ddot_vals, s2ddot)
end

# ╔═╡ 51197fff-28b2-41fa-b914-4ca2a76ce64e
md"""
	## Piston Kinematics As A Function of S1 and S2
	"""

# ╔═╡ 60f43290-b49d-43e8-981d-5d4d1e0ff2df
begin
	function piston1_xy(s1)
	    return s1*cos(α1), s1*sin(α1)
	end
	function piston2_xy(s2)
	    return -s2*cos(α2), -s2*sin(α2)
	end
end

# ╔═╡ 7230b43c-7d84-49f9-8e4c-e72302d970f7
begin
	function piston1_v_xy(s1dot)
	    return s1dot*cos(α1), s1dot*sin(α1)
	end
	function piston2_v_xy(s2dot)
	    return -s2dot*cos(α2), -s2dot*sin(α2)
	end
end

# ╔═╡ 3e362c99-b77a-49cf-a9e8-f5fc91e4d0e9
begin
	function piston1_a_xy(s1ddot)
	    return s1ddot*cos(α1), s1ddot*sin(α1)
	end
	function piston2_a_xy(s2ddot)
	    return -s2ddot*cos(α2), -s2ddot*sin(α2)
	end
end

# ╔═╡ e197d81d-799d-4671-a755-49388089028f
xp11,xp12 = piston1_xy(s1_vals)

# ╔═╡ 2987f28d-ed2a-4156-9d6c-5d9a6f98383c
xp21,xp22 = piston2_xy(s2_vals)

# ╔═╡ 8e48a543-5127-49fc-ad17-6eff85b1c0ec
begin
	xv11,xv12 = piston1_v_xy(s1dot_vals)
	xv21,xv22 = piston2_v_xy(s2dot_vals)
end

# ╔═╡ 141229c0-b611-48e2-a776-f8c539295fe8
begin
	xa11,xa12 = piston1_a_xy(s1ddot_vals)
	xa21,xa22 = piston2_a_xy(s2ddot_vals)
end

# ╔═╡ bde5d03d-0eea-4869-a93c-7c89d53b5b2d
begin
    anim = @animate for k in eachindex(tspan)

        # X position 
        p1 = plot(
            tspan, xp11,
            xlabel="Time [s]", ylabel="Position X [m]",
            title="Position X",
            lw=2, color=:steelblue, label="Piston 1"
        )
        plot!(p1, tspan, xp21, lw=2, color=:crimson, label="Piston 2")
        scatter!(p1, [tspan[k]], [xp11[k]], color=:steelblue, ms=5, label="")
        scatter!(p1, [tspan[k]], [xp21[k]], color=:crimson, ms=5, label="")
		
        # Y position
        p2 = plot(
            tspan, xp12,
            xlabel="Time [s]", ylabel="Position Y [m]",
            title="Position Y",
            lw=2, color=:steelblue, label="Piston 1"
        )
        plot!(p2, tspan, xp22, lw=2, color=:crimson, label="Piston 2")
        scatter!(p2, [tspan[k]], [xp12[k]], color=:steelblue, ms=5, label="")
        scatter!(p2, [tspan[k]], [xp22[k]], color=:crimson, ms=5, label="")

        # X velocity
      
        p3 = plot(
            tspan, xv11,
            xlabel="Time [s]", ylabel="Velocity X [m/s]",
            title="Velocity X",
            lw=2, color=:steelblue, label="Piston 1"
        )
        plot!(p3, tspan, xv21, lw=2, color=:crimson, label="Piston 2")
        scatter!(p3, [tspan[k]], [xv11[k]], color=:steelblue, ms=5, label="")
        scatter!(p3, [tspan[k]], [xv21[k]], color=:crimson, ms=5, label="")

        # Y velocity
        p4 = plot(
            tspan, xv12,
            xlabel="Time [s]", ylabel="Velocity Y [m/s]",
            title="Velocity Y",
            lw=2, color=:steelblue, label="Piston 1"
        )
        plot!(p4, tspan, xv22, lw=2, color=:crimson, label="Piston 2")
        scatter!(p4, [tspan[k]], [xv12[k]], color=:steelblue, ms=5, label="")
        scatter!(p4, [tspan[k]], [xv22[k]], color=:crimson, ms=5, label="")

        # X acceleration
        p5 = plot(
            tspan, xa11,
            xlabel="Time [s]", ylabel="Accel X [m/s²]",
            title="Acceleration X",
            lw=2, color=:steelblue, label="Piston 1"
        )
        plot!(p5, tspan, xa21, lw=2, color=:crimson, label="Piston 2")
        scatter!(p5, [tspan[k]], [xa11[k]], color=:steelblue, ms=5, label="")
        scatter!(p5, [tspan[k]], [xa21[k]], color=:crimson, ms=5, label="")

        # Y acceleration
        p6 = plot(
            tspan, xa12,
            xlabel="Time [s]", ylabel="Accel Y [m/s²]",
            title="Acceleration Y",
            lw=2, color=:steelblue, label="Piston 1"
        )
        plot!(p6, tspan, xa22, lw=2, color=:crimson, label="Piston 2")
        scatter!(p6, [tspan[k]], [xa12[k]], color=:steelblue, ms=5, label="")
        scatter!(p6, [tspan[k]], [xa22[k]], color=:crimson, ms=5, label="")

		
        plot(p1, p2, p3, p4, p5, p6; layout=(2,3), size=(1200,800))

    end

    gif(anim, "piston_kinematics_full.gif"; fps=45)
end

# ╔═╡ cb7778a3-b46b-4dac-a775-3053be621edd
begin
	 figure_p1_x = plot(
		 tspan,
		 [xp11, xv11, xa11];
		 xlabel = "Time (s)",
		 ylabel = "Vector Value",
		 label = ["Position" "Velocity" "Acceleration"],
		 linewidth = 2,
		 title = "Piston 1 X"
	 )
	
	figure_p1_y =  plot(
		 tspan,
		 [xp12, xv12, xa12];
		 xlabel = "Time (s)",
		 ylabel = "Vector Value",
		 label = ["Position" "Velocity" "Acceleration"],
		 linewidth = 2,
		 title = "Piston 1 Y"
	 )

	 figure_p2_x = plot(
		 tspan,
		 [xp21, xv21, xa21];
		 xlabel = "Time (s)",
		 ylabel = "Vector Value",
		 label = ["Position" "Velocity" "Acceleration"],
		 linewidth = 2,
		 title = "Piston 2 X"
	 )
	
	figure_p2_y =  plot(
		 tspan,
		 [xp22, xv22, xa22];
		 xlabel = "Time (s)",
		 ylabel = "Vector Value",
		 label = ["Position" "Velocity" "Acceleration"],
		 linewidth = 2,
		 title = "Piston 2 Y"
	 )


	
	plot(figure_p1_x,figure_p1_y,figure_p2_x,figure_p2_y, layout=(2,2),size=(1100,1000) ,plot_title = "Piston Kinematics")
end

# ╔═╡ e0d86b4f-cb3f-43ac-ab14-6ba1f34bffe2
begin
    Center_posx = similar(tspan)
    Center_posy = similar(tspan)

    for (i, t) in enumerate(tspan)
        Center_posx[i], Center_posy[i] = center_position(t, L)
    end
end

# ╔═╡ 05eb6373-6b67-4048-9188-05e733759b3b
begin
    # --- Bar Center X ---
    p1 = plot(
        tspan, Center_posx,
        xlabel="Time [s]", ylabel="X Position [m]",
        title="Bar Center X",
        lw=2, color=:steelblue, label="r3x"
    )

    # --- Bar Center Y ---
    p2 = plot(
        tspan, Center_posy,
        xlabel="Time [s]", ylabel="Y Position [m]",
        title="Bar Center Y",
        lw=2, color=:crimson, label="r3y"
    )

    # --- Final 2×2 Layout ---
    plot(p1, p2; layout=(1,2), size=(1200,400))
end

# ╔═╡ a1e54091-8364-47c7-8fe1-0bb52b184204
md"""# System Rotation Animation Code"""

# ╔═╡ 6dd66061-fa03-4a94-abc1-0e54152c4e93
#Defines Channel Directions for Plots
begin
	channel1_direction = (cos(α1), sin(α1))
	channel2_direction = (cos(α2), sin(α2))  # crossed channels
end

# ╔═╡ d073f553-90ea-4818-99a1-1cff4f8a5a2a
begin
	#Generates the Block Visual
	function rectangle_shape(cx, cy, dx, dy; length=0.04, width=0.02)
	    dir = [dx, dy] / norm([dx, dy])
	    nrm = [-dir[2], dir[1]]
	    hL, hW = length/2, width/2
	    pts = [ dir*hL+nrm*hW, dir*hL-nrm*hW,
	           -dir*hL-nrm*hW, -dir*hL+nrm*hW ]
	    Shape(cx .+ getindex.(pts,1), cy .+ getindex.(pts,2))
	end
	
	#Generates the slot visuals
	function slot_shape(cx, cy, dx, dy; half_length=0.25, half_width=0.025)
	    dir = [dx, dy] / norm([dx, dy])
	    nrm = [-dir[2], dir[1]]
	    ends = [-half_length*dir, half_length*dir]
	    pts = [ends[1]+half_width*nrm, ends[2]+half_width*nrm,
	           ends[2]-half_width*nrm, ends[1]-half_width*nrm]
	    Shape(cx .+ getindex.(pts,1), cy .+ getindex.(pts,2))
	end
end

# ╔═╡ 9e28d478-74a5-4ec5-ae58-5fd7890cf495
begin
	channel1 = slot_shape(0.0, 0.0, channel1_direction...; half_length=0.14)
	channel2 = slot_shape(0.0, 0.0, channel2_direction...; half_length=0.14)
end

# ╔═╡ ce209f00-f61c-4b38-9b22-fff4b5466965
begin
	time_s = tspan
	x1 = xp11; y1 = xp12
	x2 = xp21; y2 = xp22
	
	x_min = -0.15; x_max = 0.15
	y_min = -0.15; y_max = 0.15
	
	anim2 = @animate for k in eachindex(x1)
	    plot(; xlim=(x_min, x_max), ylim=(y_min, y_max), aspect_ratio=:equal,
	         legend=false, title="Mechanism motion  θ̇₃ = 2 rad/s  t=$(round(time_s[k], digits=3)) s")
	
	    plot!([x1[k], x2[k]], [y1[k], y2[k]]; lw=4)                 # connecting link
	    scatter!([x1[k]], [y1[k]]; ms=8)                             # piston 1
	    scatter!([x2[k]], [y2[k]]; ms=8)                             # piston 2
	end
	
	gif(anim2, "mechanism.gif"; fps=60)
	``
end

# ╔═╡ 5975d923-89d4-4574-a6f0-889914ab8776
dt = tspan[2]-tspan[1]

# ╔═╡ acb15deb-71aa-489e-8e15-b677b94ee0fc
begin
	anim3 = @animate for k in eachindex(x1)
	    plot(; aspect_ratio=:equal, legend=:topright,
	         xlim=(-0.35,0.35), ylim=(-0.35,0.35),
	         title="ω₃ = 2 rad/s   t=$(round(time_s[k],digits=1)) s")
	
	    plot!(channel1; c=:gray60, linecolor=:gray60, label= false)
	    plot!(channel2; c=:gray60, linecolor=:gray60, label= false)
	
	    piston1_shape = rectangle_shape(x1[k], y1[k], channel1_direction...)
	    piston2_shape = rectangle_shape(x2[k], y2[k], channel2_direction...)
	
	    plot!(piston1_shape; c=:dodgerblue, linecolor=:navy,
	          label = (k == 1 ? "Piston 1" : "Piston 1"))
	    plot!(piston2_shape; c=:tomato, linecolor=:maroon,
	          label = (k == 1 ? "Piston 2" : "Piston 2"))
	
	    plot!([x1[k], x2[k]], [y1[k], y2[k]]; lw=4, c=:black,
	          label = (k == 1 ? "Connecting bar" : "Connector Bar"))
	end
	
	gif(anim3, "mechanism_w_geometry.gif"; fps= 45)
	
end

# ╔═╡ 74b56528-a4b7-4ab8-b49f-bd1c9dd345a2
begin
	# Hanna Contribution Section
	# dependant on: L, ω3, α1, α2, tspan, s1_vals, s2_vals,
	#                s1dot_vals, s2dot_vals, s1ddot_vals, s2ddot_vals
	#                θ3 (function), xp11, xp12, xp21, xp22

	
	# contrib-01
	md"""
	---
	# 
	The sections below add four analyses not covered in the baseline notebook:
	
	1. **Constraint satisfaction check** – verifies that $C(\mathbf{q},t) = s_1^2 + s_2^2 - L^2 = 0$ holds numerically at every time step.
	2. **Speed and acceleration magnitude plots** – the baseline plots show X and Y components separately; here we show the scalar magnitudes $|\dot{\mathbf{r}}|$ and $|\ddot{\mathbf{r}}|$ for each piston.
	3. **Phase portrait** – plots $s_2$ versus $s_1$, revealing the elliptic coupling between the two pistons.
	4. **Frequency-domain analysis (FFT)** – confirms that the piston velocities are periodic at the driving frequency $\omega_3 / (2\pi)$ and identifies harmonic content.
	"""
end

# ╔═╡ da569552-04f2-47d8-8ad6-8db23c6a60c3
md"""
## 1. Constraint Satisfaction Verification

If the model is correct, $C = s_1^2 + s_2^2 - L^2$ should be identically zero (or numerically near-zero) for all $t$.
"""

# ╔═╡ 12bd10d7-df55-47f0-bd53-837ace32b56e
begin
    # Make sure all arrays have the same length
    n = length(tspan)
    constraint_residual = s1_vals[1:n].^2 .+ s2_vals[1:n].^2 .- L^2

    fig_constraint = plot(
        tspan, constraint_residual;
        xlabel = "Time (s)",
        ylabel = "C(q,t) = s₁² + s₂² - L²  [m²]",
        title  = "Constraint Residual (should be ≈ 0)",
        lw = 2,
        color  = :darkorange,
        label  = "C(q,t)",
        ylim   = (-1e-14, 1e-14),
    )
    hline!(fig_constraint, [0.0]; lw = 1, ls = :dash, color = :black, label = "zero")
    
    fig_constraint
end

# ╔═╡ fc538a2f-7fe1-4966-a945-2bf0d58f6498
begin
    using FFTW
    
    dt_val  = tspan[2] - tspan[1]
    fs      = 1.0 / dt_val
    N_fft   = n

    S1dot_fft = fft(s1dot_vals[1:n])
    freqs     = (0:N_fft-1) .* (fs / N_fft)
    amp_spec  = (2.0 / N_fft) .* abs.(S1dot_fft)

    half      = div(N_fft, 2)
    freqs_pos = freqs[1:half]
    amp_pos   = amp_spec[1:half]

    fig_fft = plot(
        freqs_pos, amp_pos;
        xlabel = "Frequency  [Hz]",
        ylabel = "Amplitude  [m/s]",
        title  = "FFT of ṡ₁  (piston 1 velocity)",
        lw = 1.5, color = :teal,
        label  = "|FFT(ṡ₁)|",
        xlim   = (0, 4.0)
    )

    f_drive = ω3 / (2π)
    vline!(fig_fft, [f_drive]; lw = 2, ls = :dash,
           color = :red, label = "f_drive ≈ $(round(f_drive, digits=3)) Hz")
    vline!(fig_fft, [2*f_drive, 3*f_drive]; lw = 1.5, ls = :dot,
           color = :darkorange, label = "harmonics")
    
    fig_fft
end

# ╔═╡ f2c7ad23-5463-4ab5-9f53-20c137cd80da
md"""
The residual sits precisely at the zero mark, confirming the analytical solution satisfies the constraint exactly. Any visible drift would indicate a numerical or modeling error.
"""

# ╔═╡ 561fcb4e-c783-4f51-beda-47247bb6eb50
md"""
## 2. Speed and Acceleration Magnitudes

The x- and y-component plots in the baseline section are useful, but the *scalar speed* $|\dot{\mathbf{r}}_i| = \sqrt{\dot{x}_i^2 + \dot{y}_i^2}$ and *scalar acceleration magnitude* give a single intuitive measure of how fast each piston moves along its track.
"""

# ╔═╡ 2b808bd4-14ef-46d7-8f81-edf5aa739db6
begin
    speed1   = abs.(s1dot_vals[1:n])
    speed2   = abs.(s2dot_vals[1:n])
    accel_mag1 = abs.(s1ddot_vals[1:n])
    accel_mag2 = abs.(s2ddot_vals[1:n])

    fig_speed = plot(
        tspan, speed1;
        xlabel = "Time (s)", ylabel = "Speed  [m/s]",
        title  = "Piston Speed Magnitudes",
        lw = 2, color = :steelblue, label = "Piston 1  |ṡ₁|"
    )
    plot!(fig_speed, tspan, speed2; lw = 2, color = :crimson, label = "Piston 2  |ṡ₂|")

    fig_accel = plot(
        tspan, accel_mag1;
        xlabel = "Time (s)", ylabel = "Acceleration  [m/s²]",
        title  = "Piston Acceleration Magnitudes",
        lw = 2, color = :steelblue, label = "Piston 1  |s̈₁|"
    )
    plot!(fig_accel, tspan, accel_mag2; lw = 2, color = :crimson, label = "Piston 2  |s̈₂|")

    plot(fig_speed, fig_accel; layout = (1, 2), size = (1100, 400),
         plot_title = "Scalar Magnitudes of Piston Kinematics")
end

# ╔═╡ 5d3d5970-2b0f-48ec-8908-d9d90a3b6527

md"""
**Observation:** both pistons reach the same peak speed and acceleration (as expected by symmetry), but they are 90° out of phase in time — when one piston is at maximum speed the other is momentarily stopped.
"""

# ╔═╡ 64138c01-10a6-4f4e-8816-d08c5beeecf5
md"""
## 3. Phase Portrait: $s_2$ vs $s_1$

A phase portrait of the two slider displacements reveals how they are geometrically coupled.
Because $s_1$ and $s_2$ are sinusoidal with the same frequency and amplitude but a 90° phase shift,
the portrait traces an ellipse (or circle in normalised units).
"""

# ╔═╡ ea59ae9f-444f-49fb-88d6-23e63c8895b6
begin
    fig_phase = plot(
        s1_vals[1:n], s2_vals[1:n];
        xlabel = "s₁  [m]",
        ylabel = "s₂  [m]",
        title  = "Phase Portrait: Piston 2 vs Piston 1 Displacement",
        lw = 2, color = :purple,
        label  = "Trajectory",
        aspect_ratio = :equal,
        legend = :topright
    )

    # Mark the starting point
    scatter!(fig_phase, [s1_vals[1]], [s2_vals[1]];
             ms = 7, color = :green, label = "t = 0")

    # Annotate quarter-period landmarks
    quarter_idx = div(n, 4)
    half_idx    = div(n, 2)
    scatter!(fig_phase,
             [s1_vals[quarter_idx], s1_vals[half_idx]],
             [s2_vals[quarter_idx], s2_vals[half_idx]];
             ms = 6, color = :darkorange,
             label = "T/4, T/2")
    
    fig_phase
end

# ╔═╡ 79a5ef9e-f0fd-497a-819e-5de25d7ab39b
md"""
The closed elliptical orbit confirms the system is conservative and periodic.
The aspect ratio of the ellipse encodes the amplitude ratio of the two pistons. The pistons  90° away from one another at all times is displayed as points on opposite ends of the closed elliptical's center.
"""

# ╔═╡ 2d4536b1-fb90-4ce3-9993-87b0a6fd4d8b
# contrib-11
md"""
## 4. Frequency-Domain Analysis

Since the bar rotates at a constant $\dot{\theta}_3 = 2$ rad/s, the fundamental driving frequency is

$$f_{\text{drive}} = \frac{\omega_3}{2\pi} \approx 0.318 \text{ Hz}$$

The piston velocities should therefore show a dominant peak at this frequency (and possibly integer harmonics).
We compute the one-sided amplitude spectrum of $\dot{s}_1$ using Julia's `FFTW` (loaded via `DSP`).
"""

# ╔═╡ 2f339107-2c86-4f8d-accc-88a5e738ec09
# contrib-fft-check
begin 
    # Find index of main frequency
    idx_main = argmin(abs.(freqs_pos .- f_drive))
    main_amp = amp_pos[idx_main]
    
    # Find max harmonic amplitude (excluding main)
    harmonic_amps = amp_pos[amp_pos .< main_amp]
    max_harmonic = maximum(harmonic_amps)
    
    ratio = max_harmonic / main_amp
    
    println("Main amplitude: $main_amp")
    println("Largest harmonic: $max_harmonic")  
    println("Harmonic ratio: $ratio")
end

# ╔═╡ adec65bf-2c38-4ec7-afed-2ae474e8cb6c
md"""
**Result:** the dominant spectral peak sits exactly at $f = \omega_3 / (2\pi) \approx 0.318$ Hz, with the next significant content at $2f$ (second harmonic). Higher harmonics decay rapidly, confirming the motion is well-described by its first two Fourier terms. This is consistent with the analytical forms $s_1(t) = -\frac{L}{2\sqrt{2}}(\cos\theta_3 + \sin\theta_3)$ which contains a single sinusoidal frequency. Harmonic content is <1% of the fundamental, confirming near-perfect sinusoidal motion.
"""

# ╔═╡ c167d7fe-18a2-490d-bf61-492cef7cda38
md"""
## Summary

| Analysis | Key findings |
|---|---|
| Constraint residual | analytical solution satisfies constraint (residiual ≈ 0) |
| Speed magnitudes | Pistons reach equal peak speeds, 90° out of phase |
| Phase portrait | Closed ellipse confirms conservative, periodic motion |
| FFT | Single dominant frequency at $\omega_3/(2\pi)$; higher harmonics negligible |
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
DSP = "717857b8-e6f2-59f4-9121-6e50c889abd2"
FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
InteractiveUtils = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Markdown = "d6f4376e-aef5-505a-96c1-9c027394607a"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"

[compat]
DSP = "~0.8.4"
FFTW = "~1.10.0"
Plots = "~1.41.6"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.10"
manifest_format = "2.0"
project_hash = "4db9c35951d8b12f93612f55ccc633d948a79d5b"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

    [deps.AbstractFFTs.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d0efe2c6fdcdaa1c161d206aa8b933788397ec71"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.6+0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b0fd3f56fa442f81e0a47815c92245acfaaa4e34"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.31.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "67e11ee83a43eb71ddc950302c53bf33f0690dfe"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.12.1"

    [deps.ColorTypes.extensions]
    StyledStringsExt = "StyledStrings"

    [deps.ColorTypes.weakdeps]
    StyledStrings = "f489334b-da3d-4c2e-b8f0-e476e12c162b"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "8b3b6f87ce8f65a2b4f857528fd8d70086cd72b1"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.11.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "37ea44092930b1811e666c3bc38065d7d87fcc74"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.13.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "21d088c496ea22914fe80906eb5bce65755e5ec8"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.1"

[[deps.ConstructionBase]]
git-tree-sha1 = "b4b092499347b18a015186eae3042f72267106cb"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.6.0"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DSP]]
deps = ["Bessels", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "5989debfc3b38f736e69724818210c67ffee4352"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.8.4"

    [deps.DSP.extensions]
    OffsetArraysExt = "OffsetArrays"

    [deps.DSP.weakdeps]
    OffsetArrays = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["OrderedCollections"]
git-tree-sha1 = "e86f4a2805f7f19bec5129bc9150c38208e5dc23"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.19.4"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "473e9afc9cf30814eb67ffa5f2db7df82c3ad9fd"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.16.2+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
git-tree-sha1 = "7442a5dfe1ebb773c29cc2962a8980f47221d76c"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.5"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "27af30de8b5445644e8ffe3bcb0d72049c089cf1"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.7.3+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "95ecf07c2eea562b5adbd0696af6db62c0f52560"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.5"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libva_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "66381d7059b5f3f6162f28831854008040a4e905"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "8.0.1+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "Libdl", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "97f08406df914023af55ade2f843c39e99c5d969"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.10.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6d6219a004b8cf1e0b4dbe27a2860b8e04eba0be"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.11+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "f85dac9a96a01087df6e3a749840015a0ca3817d"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.17.1+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "70329abc09b886fd2c5d94ad2d9527639c421e3e"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.14.3+1"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7a214fdac5ed5f59a22c2d9a885a16da1c74bbc7"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.17+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "b7bfd56fa66616138dfe5237da4dc13bbd83c67f"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.1+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "44716a1a667cb867ee0e9ec8edc31c3e4aa5afdc"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.24"

    [deps.GR.extensions]
    IJuliaExt = "IJulia"

    [deps.GR.weakdeps]
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "be8a1b8065959e24fdc1b51402f39f3b6f0f6653"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.24+0"

[[deps.GettextRuntime_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll"]
git-tree-sha1 = "45288942190db7c5f760f59c04495064eedf9340"
uuid = "b0724c58-0f36-5564-988d-3bb0596ebc4a"
version = "0.22.4+0"

[[deps.Ghostscript_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Zlib_jll"]
git-tree-sha1 = "38044a04637976140074d0b0621c1edf0eb531fd"
uuid = "61579ee1-b43e-5ca0-a5da-69d92c66a64b"
version = "9.55.1+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "GettextRuntime_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "24f6def62397474a297bfcec22384101609142ed"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.86.3+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "51059d23c8bb67911a2e6fd5130229113735fc7e"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.11.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "f923f9a774fcf3f5cb761bfa43aeadd689714813"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.1+0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "ec1debd61c300961f98064cfb21287613ad7f303"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "b2d91fe939cae05960e760110b328288867b5758"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.6"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "82f7acdc599b65e0f8ccd270ffa1467c21cb647b"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.11"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "0533e564aae234aff59ab625543145446d8b6ec2"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.1"

[[deps.JSON]]
deps = ["Dates", "Logging", "Parsers", "PrecompileTools", "StructUtils", "UUIDs", "Unicode"]
git-tree-sha1 = "67c6f1f085cb2671c93fe34244c9cccde30f7a26"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "1.5.0"

    [deps.JSON.extensions]
    JSONArrowExt = ["ArrowTypes"]

    [deps.JSON.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6893345fd6658c8e475d40155789f4860ac3b21"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.4+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "059aabebaa7c82ccb853dd4a0ee9d17796f7e1bc"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.3+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eb62a3deb62fc6d8822c0c4bef73e4412419c5d8"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.8+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "Ghostscript_jll", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "44f93c47f9cd6c7e431f2f2091fcba8f01cd7e8f"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.10"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"
    TectonicExt = "tectonic_jll"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"
    tectonic_jll = "d7dd28d6-a5e6-559c-9131-7eb760cdacc5"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c8da7e6a91781c41a863611c7e966098d783c57a"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.4.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "d36c21b9e7c172a44a10484125024495e2625ac0"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.1+1"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "97bbca976196f2a1eb9607131cb108c69ec3f8a6"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "f04133fe05eff1667d2054c53d59f9122383fe05"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.2+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d0205286d9eceadc518742860bf23f703779a3d6"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.3+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f00544d95982ea270145636c181ceda21c4e2575"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.2.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "282cadc186e7b2ae0eeadbd7a4dffed4196ae2aa"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.2.0+0"

[[deps.MacroTools]]
git-tree-sha1 = "1e0228a030642014fe5cfe68c2c0a818f9e3f522"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.16"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "8785729fa736197687541f7053f6d8ab7fc44f92"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.10"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "b513cedd20d9c914783d8ad83d08120702bf2c77"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.3"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "9b8215b1ee9e78a293f99797cd31375471b2bcae"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.3"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6aa4566bb7ae78498a5e68943863fa8b5231b59"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.6+0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.5+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "NetworkOptions", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "1d1aaa7d449b58415f97d2839c318b70ffb525a0"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.6.1"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "2ac022577e5eac7da040de17776d51bb770cd895"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.5.6+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e2bb57a313a74b8104064b7efd01406c0a50d2ff"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.6.1+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "05868e21324cede2207c6f0f466b4bfef6d5e7ee"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.1"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0662b083e11420952f2e62e17eddae7fc07d5997"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.57.0+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "7d2f8f21da5db6a806faf7b9b292296da42b2810"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.3"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "26ca162858917496748aad52bb5d3be4d26a228a"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.4"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "cb20a4eacda080e517e4deb9cfb6c7c518131265"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.41.6"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "Setfield", "SparseArrays"]
git-tree-sha1 = "2d99b4c8a7845ab1342921733fa29366dae28b24"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.1.1"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieExt = "Makie"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"
    PolynomialsRecipesBaseExt = "RecipesBase"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
    RecipesBase = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "8b770b60760d4451834fe79dd483e318eee709c4"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PtrArrays]]
git-tree-sha1 = "4fbbafbc6251b883f4d2705356f3641f3652a7fe"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.4.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "d7a4bff94f42208ce3cf6bc8e4e7d1d663e7ee8b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.10.2+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll", "Qt6Svg_jll"]
git-tree-sha1 = "d5b7dd0e226774cbd87e2790e34def09245c7eab"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.10.2+1"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "4d85eedf69d875982c46643f6b4f66919d7e157b"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.10.2+1"

[[deps.Qt6Svg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "81587ff5ff25a4e1115ce191e36285ede0334c9d"
uuid = "6de9746b-f93d-5813-b365-ba18ad4a9cf3"
version = "6.10.2+0"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "672c938b4b4e3e0169a07a5f227029d4905456f2"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.10.2+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "9b81b8393e50b7d4e6d0a9f14e192294d3b7c109"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.3.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "64d974c2e6fdf07f8155b5b2ca2ffa9069b608d9"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.2"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2700b235561b0335d5bef7097a111dc513b8655e"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.7.2"

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

    [deps.SpecialFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "4f96c596b8c8258cc7d3b19797854d368f243ddc"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.4"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6ab403037779dae8c514bad259f32a447262455a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.4"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "178ed29fd5b2a2cfc3bd31c13375ae925623ff36"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.8.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "IrrationalConstants", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "aceda6f4e598d331548e04cc6b2124a6148138e3"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.10"

[[deps.StructUtils]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "fa95b3b097bcef5845c142ea2e085f1b2591e92c"
uuid = "ec057cc2-7a8d-4b58-b3b3-92acb9f63b42"
version = "2.7.1"

    [deps.StructUtils.extensions]
    StructUtilsMeasurementsExt = ["Measurements"]
    StructUtilsStaticArraysCoreExt = ["StaticArraysCore"]
    StructUtilsTablesExt = ["Tables"]

    [deps.StructUtils.weakdeps]
    Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tables = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.URIs]]
git-tree-sha1 = "bef26fb046d031353ef97a82e3fdb6afe7f21b1a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.6.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "96478df35bbc2f3e1e791bc7a3d0eeee559e60e9"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.24.0+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b29c22e245d092b8b4e8d3c09ad7baa586d9f573"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.8.3+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a3ea76ee3f4facd7a64684f9af25310825ee3668"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.2+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "9c7ad99c629a44f81e7799eb05ec2746abb5d588"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.6+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "808090ede1d41644447dd5cbafced4731c56bd2f"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.13+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aa1261ebbac3ccc8d16558ae6799524c450ed16b"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.13+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c74ca84bbabc18c4547014765d194ff0b4dc9da"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.4+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "52858d64353db33a56e13c341d7bf44cd0d7b309"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.6+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "1a4a26870bf1e5d26cd585e38038d399d7e65706"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.8+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "75e00946e43621e09d431d9b95818ee751e6b2ef"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.2+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "a376af5c7ae60d29825164db40787f15c80c7c54"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.3+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "0ba01bc7396896a4ace8aab67db31403c71628f4"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.7+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "6c174ef70c96c76f4c3f4d3cfbe09d018bcd1b53"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.6+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "7ed9347888fac59a618302ee38216dd0379c480d"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.12+0"

[[deps.Xorg_libpciaccess_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "4909eb8f1cbf6bd4b1c30dd18b2ead9019ef2fad"
uuid = "a65dc6b1-eb27-53a1-bb3e-dea574b5389e"
version = "0.18.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXau_jll", "Xorg_libXdmcp_jll"]
git-tree-sha1 = "bfcaf7ec088eaba362093393fe11aa141fa15422"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.1+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "ed756a03e95fff88d8f738ebc2849431bdd4fd1a"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.2.0+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "9750dc53819eba4e9a20be42349a6d3b86c7cdf8"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.6+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f4fc02e384b74418679983a97385644b67e1263b"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll"]
git-tree-sha1 = "68da27247e7d8d8dafd1fcf0c3654ad6506f5f97"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "44ec54b0e2acd408b0fb361e1e9244c60c9c3dd4"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.1+0"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "5b0263b6d080716a02544c55fdff2c8d7f9a16a0"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.10+0"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_jll"]
git-tree-sha1 = "f233c83cad1fa0e70b7771e0e21b061a116f2763"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.2+0"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "801a858fc9fb90c11ffddee1801bb06a738bda9b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.7+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "00af7ebdc563c9217ecc67776d1bbf037dbcebf4"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.44.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c3b0e6196d50eab0c5ed34021aaa0bb463489510"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.14+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b6a34e0e0960190ac2a4363a1bd003504772d631"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.61.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "371cc681c00a3ccc3fbc5c0fb91f58ba9bec1ecf"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.13.1+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "125eedcb0a4a0bba65b657251ce1d27c8714e9d6"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.17.4+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libdrm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libpciaccess_jll"]
git-tree-sha1 = "63aac0bcb0b582e11bad965cef4a689905456c03"
uuid = "8e53e030-5e6c-5a89-a30b-be5b7263a166"
version = "2.4.125+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56d643b57b188d30cccc25e331d416d3d358e557"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.13.4+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "646634dd19587a56ee2f1199563ec056c5f228df"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.4+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "91d05d7f4a9f67205bd6cf395e488009fe85b499"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.28.1+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e2a7072fc0cdd7949528c1455a3e5da4122e1153"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.56+0"

[[deps.libva_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll", "Xorg_libXfixes_jll", "libdrm_jll"]
git-tree-sha1 = "7dbf96baae3310fe2fa0df0ccbb3c6288d5816c9"
uuid = "9a156e7d-b971-5f62-b2c9-67348b8fb97c"
version = "2.23.0+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll"]
git-tree-sha1 = "11e1772e7f3cc987e9d3de991dd4f6b2602663a5"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.8+0"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "b4d631fd51f2e9cdd93724ae25b2efc198b059b1"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.7+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "1350188a69a6e46f799d3945beef36435ed7262f"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14cc7083fc6dff3cc44f2bc435ee96d06ed79aa7"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "10164.0.1+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e7b67590c14d487e734dcb925924c5dc43ec85f3"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "4.1.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "a1fc6507a40bf504527d0d4067d718f8e179b2b8"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.13.0+0"
"""

# ╔═╡ Cell order:
# ╠═f5b9dd00-5682-4241-8b80-62ab0ac3bda8
# ╠═5f84fe8f-ee3c-4963-b458-5d5b11450f5d
# ╠═6ce190b9-7344-4cbf-8e12-adefc8a3bba2
# ╠═943023c9-f07b-486d-9b0e-3983b0147af6
# ╠═758462d7-cf01-4f41-9641-b806a3adfbfe
# ╠═c46aa017-b129-4596-aa65-15d4d4aca963
# ╠═8cba0c20-1450-47ce-b2eb-714f18916f4f
# ╠═1b810674-e02c-452e-b32a-d4137b39ea8e
# ╠═694a9ed8-8914-4571-8298-505dd11f6b78
# ╠═370a9b3c-236d-46c7-8d6f-4ee25947ec61
# ╠═1e7c0d64-3b10-4162-86c2-9c6b15ec644d
# ╠═73a6416f-506f-4eec-a03f-a7a9bec789d0
# ╠═666b9581-f2f9-4c58-b245-59ad60803772
# ╠═179b675e-cad5-47c2-9df1-0c6db8541b79
# ╠═ef347f87-3dbc-4fbd-bef5-8008daccfa7a
# ╠═a69048f1-fc46-412a-8455-7e7b8048c7ba
# ╠═0e457baa-f346-427e-a5dc-6f2c8be711eb
# ╠═72752148-5b7d-454e-b377-0cfaae37b33f
# ╠═e8799a11-75a5-4446-9518-d5f2f909c283
# ╠═f8e5cfa9-3765-43e0-9bae-dc245f7c1eff
# ╠═5d2dcd83-64e9-40b1-8e3e-3d3150d06e76
# ╠═dcd74d4d-44a8-4ab7-815f-aa7bed0da8d2
# ╠═a8e64e00-b8d2-4cc6-b2c9-148e84364490
# ╠═e3a09adf-a08a-4e38-8245-87d533038b13
# ╠═e2836679-2a9e-48fa-8a1b-dd658d07587a
# ╠═51197fff-28b2-41fa-b914-4ca2a76ce64e
# ╠═60f43290-b49d-43e8-981d-5d4d1e0ff2df
# ╠═7230b43c-7d84-49f9-8e4c-e72302d970f7
# ╠═3e362c99-b77a-49cf-a9e8-f5fc91e4d0e9
# ╠═e197d81d-799d-4671-a755-49388089028f
# ╠═2987f28d-ed2a-4156-9d6c-5d9a6f98383c
# ╠═8e48a543-5127-49fc-ad17-6eff85b1c0ec
# ╠═141229c0-b611-48e2-a776-f8c539295fe8
# ╠═bde5d03d-0eea-4869-a93c-7c89d53b5b2d
# ╠═cb7778a3-b46b-4dac-a775-3053be621edd
# ╠═e0d86b4f-cb3f-43ac-ab14-6ba1f34bffe2
# ╠═05eb6373-6b67-4048-9188-05e733759b3b
# ╠═a1e54091-8364-47c7-8fe1-0bb52b184204
# ╠═6dd66061-fa03-4a94-abc1-0e54152c4e93
# ╠═d073f553-90ea-4818-99a1-1cff4f8a5a2a
# ╠═9e28d478-74a5-4ec5-ae58-5fd7890cf495
# ╠═ce209f00-f61c-4b38-9b22-fff4b5466965
# ╠═5975d923-89d4-4574-a6f0-889914ab8776
# ╠═acb15deb-71aa-489e-8e15-b677b94ee0fc
# ╠═74b56528-a4b7-4ab8-b49f-bd1c9dd345a2
# ╠═da569552-04f2-47d8-8ad6-8db23c6a60c3
# ╠═12bd10d7-df55-47f0-bd53-837ace32b56e
# ╠═f2c7ad23-5463-4ab5-9f53-20c137cd80da
# ╠═561fcb4e-c783-4f51-beda-47247bb6eb50
# ╠═2b808bd4-14ef-46d7-8f81-edf5aa739db6
# ╠═5d3d5970-2b0f-48ec-8908-d9d90a3b6527
# ╠═64138c01-10a6-4f4e-8816-d08c5beeecf5
# ╠═ea59ae9f-444f-49fb-88d6-23e63c8895b6
# ╠═79a5ef9e-f0fd-497a-819e-5de25d7ab39b
# ╠═2d4536b1-fb90-4ce3-9993-87b0a6fd4d8b
# ╠═fc538a2f-7fe1-4966-a945-2bf0d58f6498
# ╠═2f339107-2c86-4f8d-accc-88a5e738ec09
# ╠═adec65bf-2c38-4ec7-afed-2ae474e8cb6c
# ╠═c167d7fe-18a2-490d-bf61-492cef7cda38
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
