import numpy as np
from graphviz import Digraph

DIMENSION=1

smoothing = 0.001
G = 6.67430e-4	 # gravitational constant in m^3 kg^-1 s^-2

def grav(p1, p2, m1, m2):
	r = p1 - p2
	distance = np.sqrt(np.linalg.norm(r)**2 + smoothing**2)
	return G * m1 * m2 * r / distance**3
				

class TreeNode:
	def __init__(self, value, center, length):
		self.state_vector = value
		self.center = center
		self.length = length
		self.tm = 0.0
		self.cm = 0.0
		self.left = None
		self.right = None

class BinaryTree:
	depth = 0
	def __init__(self, box_max):		
		self.root = TreeNode(None, 0., box_max)
		self.root.left = TreeNode(None, -box_max/2, box_max/2)
		self.root.right = TreeNode(None, box_max/2, box_max/2)

	def _internal(self, node):
		return node.left is not None or node.right is not None
	def _external(self, node):
		return node.state_vector is not None and node.left is None and node.right is None
	def _empty_value(self, node):
		return node.state_vector is None
	
	"""
	If node x is an internal node, update the center-of-mass and total mass of x.
	Recursively insert the body b in the appropriate quadrant.
	"""
	def _internal_insert(self, node: TreeNode, value: float):
		if node.center > value[X_NDX]:
			self._insert(node.left, value)
		else:
			self._insert(node.right, value)
		
			
	"""
	If node x is an external node, say containing a body named c,
	then there are two bodies b and c in the same region. Subdivide the region
	further by creating four children. Then, recursively insert both b and c into the
	appropriate quadrant(s). Since b and c may still end up in the same quadrant,
	there may be several subdivisions during a single insertion. Finally, update the
	center-of-mass and total mass of x.
	"""

	def _update_node(self, node, value):
		#node.cm = (node.cm*node.tm + value.mass*value.position)/(node.tm + value.mass)
		node.tm += value.mass
		
	
	def _external_insert(self, node, value):
		old_val = node.state_vector
		node.value = None
		node.left = TreeNode(None, node.center - node.length/2, node.length/2)
		node.right = TreeNode(None, node.center + node.length/2, node.length/2)
		self.insert(value)
		self.insert(old_val)
		
	def insert(self, value):
		if value[X_NDX] > self.root.length + self.root.center or value[X_NDX] < -1*self.root.length - self.root.center:
			print('Invalid value in insertion')
			exit(-1)
		
		self._insert(self.root, value)
	
	def _insert(self, node, value):
		if node is None:
			print('uh oh')
			return
		
		if self._internal(node):
			self._internal_insert(node, value)
		elif self._external(node):
			self._external_insert(node, value)
		elif self._empty_value(node):
			node.state_vector = value
			#self._update_node(node, value)
		else:
			print('Uncaught case :(')
			exit(-1)
	"""
	If the current node is an external node (and it is not body b),
	calculate the force exerted by the current node on b, and add this amount to b’s net force.
	
	Otherwise, calculate the ratio s/d. If s/d < θ, treat this internal node as a single body,
	and calculate the force it exerts on body b, and add this amount to b’s net force.

	Otherwise, run the procedure recursively on each of the current node’s children.
	"""

	def _walk_force(self, node, value, theta):
		v_position = value[X_NDX]
		v_mass = value[M_NDX]

		if node is None:
			return

		if self._external(node) and not state_equal(value, node.state_vector):
			n_position = node.state_vector[X_NDX]
			n_mass = node.state_vector[M_NDX]
			value[F_NDX] -= grav(v_position, n_position,
									 v_mass, n_mass)
			return

		dist = abs(node.cm-value[X_NDX])
		if self._internal(node) and not dist == 0 and node.length/dist < theta:
			n_position = node.cm
			n_mass = node.tm
			value[F_NDX] -= grav(v_position, n_position,
									 v_mass, n_mass)
			return
		else:
			self._walk_force(node.left, value, theta)
			self._walk_force(node.right, value, theta)
		
	def walk_force(self, value, theta):
		self._walk_force(self.root, value, theta)
	
	def _update_tm(self, node):
		if self._external(node) and node.state_vector is not None:
			return (node.state_vector[M_NDX], node.state_vector[X_NDX])
		elif self._internal(node):
			left = self._update_tm(node.left)
			right = self._update_tm(node.right)
			node.cm = (left[0]*left[1] + right[0]*right[1])/(left[0] + right[0])
			node.tm += left[0]
			node.tm += right[0]
			return (node.tm, node.cm)
		elif self._empty_value(node):
			return (0, 0)
		else:
			return (0, 0)
	def update_tm(self):
		self._update_tm(self.root)

		
		
		"""
		if self.root.value is None: #Case 1 empty node
			self.root = TreeNode(value, value, 10000)
		else:
			self._insert_recursive(self.root, value)
		"""
					
	def inorder_traversal(self, node, result=None):
		if result is None:
			result = []
		if node:
			self.inorder_traversal(node.left, result)
			if not node.state_vector is None:
				result.append(node.state_vector)
			self.inorder_traversal(node.right, result)
		return result

	def visualize_tree(self):
		if self.root is None:
			return 'The tree is empty'
		dot = Digraph()
		self._add_nodes_edges(self.root, dot)
		return dot

	def _node_string(self, node):
		return str(node)
	
	def _add_nodes_edges(self, node, dot):
		if node:
			root = self._node_string(node)
			dot.node(root)
			if node.left:
				left = self._node_string(node.left)
				dot.node(left)
				dot.edge(root, left)
				self._add_nodes_edges(node.left, dot)
			if node.right:
				right = self._node_string(node.right)
				dot.node(right)
				dot.edge(root, right)
				self._add_nodes_edges(node.right, dot)
			

# Create a binary tree instance and insert elements
	
def create_tree(positions):
	tree = BinaryTree(500)
	masses = [100000]*len(positions)
	for p, m in zip(positions, masses):
		
		tree.insert(np.array([p, 0, 0, m]))
	tree.update_tm()
	return tree

def iterate():
	nodes = T.inorder_traversal(T.root)
	for node in nodes:
		#Reset Acceleration to 0
		node.acceleration = 0
		T.walk_force(node, 0.5)
		node.acceleration /= node.mass
		node.velocity += node.acceleration * dt


	
def calculate_kinetic_energy(bodies):
	total_ke = 0.0
	for body in bodies:
		total_ke += 0.5 * body[-1] * body[1]**2
	return total_ke

def calculate_potential_energy(bodies):
	total_pe = 0.0
	n = len(bodies)
	for i in range(n):
		for j in range(i + 1, n):
			r = abs(bodies[i][0] - bodies[j][0])
			total_pe -= G * bodies[i][-1] * bodies[j][-1] / r
	return total_pe

def calculate_total_energy(bodies):
	ke = calculate_kinetic_energy(bodies)
	pe = calculate_potential_energy(bodies)
	return ke + pe


def main_loop(n, theta):
	N = len(Bodies)
	BoxSize = 2*15000*100
	
	Ke = calculate_kinetic_energy(Bodies)
	Ge = calculate_potential_energy(Bodies)
	initial_energy = Ke+Ge
	for iter_ndx in range(n):
		max_val = 0.0
		Tree = BinaryTree(BoxSize)
		for b in Bodies:
			#print(b.position)
			Tree.insert(b)
			Tree.update_tm()
		
		for b in Bodies:
			Tree.walk_force(b, theta)
		for nn in range(N):
			Bodies[nn][0] += Bodies[nn][1]*dt + 0.5*Bodies[nn][2]/Bodies[nn][3]*dt*dt
			Bodies[nn][1] += Bodies[nn][2]*dt/Bodies[nn][3]
			Bodies[nn][2] = 0.0

			positions[iter_ndx, nn] = Bodies[nn][0]
			
			if abs(Bodies[nn][0]) > max_val:
				max_val = abs(Bodies[nn][0])		
		BoxSize = 2*max_val
		"""
		if iter_ndx % 1000 == 0:
			print('ndx: ', iter_ndx)
			Ke = calculate_kinetic_energy(Bodies)
			Ge = calculate_potential_energy(Bodies)
			print(Ke+Ge, initial_energy)
		"""
X_NDX=0
V_NDX=1
F_NDX=2
M_NDX=3
#Position velocity force mass
def PointMass(position, mass):
	return np.array([position, 0., 0., mass], dtype=np.float64)


def state_equal(s1, s2):
	return s1[0] == s2[0] and s1[1] == s2[1] and s1[2] == s2[2] and s1[3] == s2[3]


Bodies = [PointMass(10., 1.), PointMass(-10., 1.), PointMass(-30., 1.), PointMass(40., 1.0)]
#Bodies = [0]*100
#for ndx in range(0, 100):
#	Bodies[ndx] = PointMass(ndx*15000, 100e3)
	
dt = 0.01
STEPS = 10000
positions = np.zeros((STEPS, len(Bodies))) 
main_loop(STEPS, 0.0)
np.save('meow.npy', positions)
