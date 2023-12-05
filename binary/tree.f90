program tree_sim
  implicit none

  type TreeNode
     real(8) :: Center  ! Center of the line segment
     real(8) :: length           ! Length of the line segment
     integer :: state_vector_index !Bootleg pointer lmfao.
     real(8) :: total_mass
     real(8) :: center_mass     
     type (TreeNode), pointer :: left, right
  end type TreeNode
  
  type BinaryTree
     type (TreeNode), pointer :: root     
  end type BinaryTree

  type (BinaryTree) :: Tree
  integer, parameter :: system_size=4
  real(8), parameter :: theta = 0.0
  real(8), parameter :: eps = 0.001
  real(8), parameter :: G =  6.67430d-4 !gravitational constant in m^3 kg^-1 s^-2
  real(8), parameter :: dt = 0.01
  
  real(8) :: STATES(system_size, 4)  ! (X, Vx, Fx, Mass)
  integer :: n, nn

  STATES(1, :) = (/  10.0d0, 0.0d0, 0.0d0, 1.0d0 /)
  STATES(2, :) = (/ -10.0d0, 0.0d0, 0.0d0, 1.0d0 /)
  STATES(3, :) = (/ -30.0d0, 0.0d0, 0.0d0, 1.0d0 /)
  STATES(4, :) = (/  40.0d0, 0.0d0, 0.0d0, 1.0d0 /)
  open(file='tree_test.dat', unit=10)
  do nn = 1, 10000
     call CREATE_TREE(Tree)
     call UPDATE_VIRTUAL_NODES(Tree)
     do n=1, system_size
        STATES(n, 3) = calculate_force(Tree, n) 
     enddo
     STATES(:, 1) = STATES(:, 1) + STATES(:, 2)*dt + 0.5*STATES(:, 3)/STATES(:, 4) *dt *dt
     STATES(:, 2) = STATES(:, 2) +  dt * STATES(:, 3)/STATES(:, 4)
     STATES(:, 3) = 0.0
     write(10,*) STATES(:, 1)
     !call WRITE_TREE(Tree % root)
     call DESTROY_TREE(Tree)
  enddo
  close(unit=10)
contains


  real(8) function grav(p1, p2, m1, m2)
    real(8) :: p1, p2, m1, m2
    real(8) :: r, distance
    r = p1 - p2
    distance = sqrt(r**2 + eps**2)
    grav = G * m1 * m2 * r / distance**3
  end function grav
    
  
  recursive real(8) function calculate_force_(Tree, node, state_vector_index) RESULT(force)
    type (BinaryTree) :: Tree
    type (TreeNode), pointer :: node
    integer :: state_vector_index
    real(8) :: dist

    if (.not. associated(node)) then
       force = 0.0
    else
         
       dist = abs(node % center_mass - STATES(state_vector_index, 1))
       if (external_node(node) .and. .not. state_vector_index == node % state_vector_index) then
          force = -1.0*grav(STATES(state_vector_index, 1),  &
               STATES(node % state_vector_index, 1), &
               STATES(state_vector_index, 4), &
               STATES(node % state_vector_index, 4))
       else if(internal_node(node) .and. dist > 0.1 .and. node % length/dist < theta) then
          force = -1.0*grav(STATES(state_vector_index, 1), &
               node % center_mass, STATES(state_vector_index, 4), node % total_mass)
       else
          force = calculate_force_(Tree, node % left, state_vector_index) + &
               calculate_force_(Tree, node % right, state_vector_index)
       endif
    endif
  end function calculate_force_
  
  recursive  real(8) function calculate_force(Tree, state_vector_index)
    type (BinaryTree) :: Tree
    integer :: state_vector_index
    calculate_force = CALCULATE_FORCE_(Tree, Tree % root, state_vector_index)    
  end function calculate_force

  SUBROUTINE CREATE_TREE(Tree)
    type (BinaryTree) :: Tree
    integer :: n
    call INIT_TREE(Tree, 0.0d0, maxval((/ maxval(STATES(:, 1)), abs(minval(STATES(:, 1))) /)))
    do n=1, system_size
       call INSERT(Tree, n)  
    enddo
  end SUBROUTINE CREATE_TREE
  
  RECURSIVE SUBROUTINE DESTROY_TREE_(Tree, Node)
    type (BinaryTree) :: Tree
    type (TreeNode), pointer :: node
    if (associated(node)) then
       call DESTROY_TREE_(Tree, node % left)
       call DESTROY_TREE_(Tree, node % right)
       deallocate(node)       
    endif       
  END SUBROUTINE DESTROY_TREE_
  
  RECURSIVE SUBROUTINE DESTROY_TREE(Tree)
    type (BinaryTree) :: Tree
    call DESTROY_TREE_(Tree, Tree % root)
  END SUBROUTINE DESTROY_TREE
 
    
  RECURSIVE SUBROUTINE UPDATE_VIRTUAL_NODES_(Tree, node)
    type (BinaryTree) :: Tree
    type (TreeNode) :: node
    if (external_node(node) .and. node % state_vector_index > 0) then

       node % total_mass = STATES(node % state_vector_index, 4)
       node % center_mass = STATES(node % state_vector_index, 1)       
    else if (internal_node(node)) then
       call UPDATE_VIRTUAL_NODES_(Tree, node % left)
       call UPDATE_VIRTUAL_NODES_(Tree, node % right)
       node % center_mass = (node % left % center_mass * node % left % total_mass + &
            node % right % center_mass * node % right % total_mass)/ &
            (node % left % total_mass + node % right % total_mass)
       node % total_mass = node % left % total_mass + node % right % total_mass
    endif
  END SUBROUTINE UPDATE_VIRTUAL_NODES_
  
  RECURSIVE SUBROUTINE UPDATE_VIRTUAL_NODES(Tree)
    type (BinaryTree) :: Tree
    call UPDATE_VIRTUAL_NODES_(Tree, Tree % root)
  END SUBROUTINE UPDATE_VIRTUAL_NODES

  SUBROUTINE WRITE_NODE(node)
    type (TreeNode), pointer :: node
    write(*,*) 'center', node % center, 'length', node % length, &
         'total mass', node % total_mass, 'center mass', &
         node % center_mass, 'LEAF', node % state_vector_index > 0 
  END SUBROUTINE WRITE_NODE
     
  RECURSIVE SUBROUTINE WRITE_TREE(node)
    type (TreeNode), pointer :: node
    if (associated(node)) then
       call WRITE_TREE(node % left)
       call WRITE_NODE(node)
       call WRITE_TREE(node % right)      
    endif    
  end SUBROUTINE WRITE_TREE
  
  RECURSIVE SUBROUTINE WRITE_LEAVES(node)
    type (TreeNode), pointer :: node
    if (associated(node)) then
       call WRITE_TREE(node % left)
       if (node % state_vector_index > 0) then
          write(*,*) STATES(node % state_vector_index, :)
       endif
       call WRITE_TREE(node % right)      
    endif    
  end SUBROUTINE WRITE_LEAVES


  logical function internal_node(Node)
    type (TreeNode) :: Node
    internal_node = associated(Node % left) .or. associated(Node % right)
  end function internal_node

  
  logical function external_node(Node)
    type (TreeNode) :: Node
    external_node = .not. associated(Node % left) .and. .not. associated(Node % right)
  end function external_node

  RECURSIVE subroutine INTERNAL_INSERT_(Tree, current, state_vector_index)
    type (BinaryTree) :: Tree
    type (TreeNode) :: current
    integer :: state_vector_index

    if (current % center > STATES(state_vector_index, 1)) then
       call INSERT_(Tree, current % left, state_vector_index)
    else
       call INSERT_(Tree, current % right, state_vector_index)
    endif    
  end subroutine INTERNAL_INSERT_

  RECURSIVE subroutine EXTERNAL_INSERT_(Tree, current, state_vector_index)
    type (BinaryTree) :: Tree
    type (TreeNode) :: current
    integer :: state_vector_index
    integer :: old
    real(8) :: center, length
    
    old = current % state_vector_index
    current % state_vector_index = -1 !(/ 0.0, 0.0, 0.0, -1.0 /)
    center = current % center
    length = current % length
    
    call Create_Node(current % left, center-length/2., length/2.)
    call Create_Node(current % right, center+length/2., length/2.) 

    call INSERT_(Tree, Tree % root, state_vector_index)
    call INSERT_(Tree, Tree % root, old)    
  end subroutine EXTERNAL_INSERT_

  RECURSIVE subroutine INSERT_(Tree, current, state_vector_index)
    type (BinaryTree) :: Tree
    type (TreeNode) :: current
    integer :: state_vector_index

    if (internal_node(current)) then
       call INTERNAL_INSERT_(Tree, current, state_vector_index)
    else if(external_node(current)) then
       if (current % state_vector_index < 0) then
          current % state_vector_index = state_vector_index
       else
          call EXTERNAL_INSERT_(Tree, current, state_vector_index)
       endif
    else
       write(*,*) 'Uncaught case in insert'
    endif
  end subroutine INSERT_

  subroutine INSERT(Tree, state_vector_index)
    type (BinaryTree) :: Tree
    integer :: state_vector_index
    real(8) :: bound_left, bound_right

    bound_left = Tree % root % center - Tree % root % length
    bound_right = Tree % root % center + Tree % root % length
    
    
    if (STATES(state_vector_index, 1) < bound_left .or. &
         STATES(state_vector_index, 1) > bound_right) then
       write(*,*) 'point outside of tree bounds... Exiting.'
       call exit(-1)
    else
       call INSERT_(Tree, Tree % root, state_vector_index)
    endif
    
  end subroutine INSERT
  
  subroutine CREATE_NODE(node, center, length)
    type (TreeNode), pointer :: node
    real(8) :: center, length
    allocate (node)
    node % center = center
    node % length = length
    node % state_vector_index = -1
    nullify (node % left)
    nullify (node % right)    
  end subroutine CREATE_NODE
  
  subroutine INIT_TREE(Tree, center, length)
    real(8) :: center, length
    type (BinaryTree) :: Tree
    call Create_Node(Tree % root, center, length)
    call Create_Node(Tree % root % left, center-length/2., length/2.)
    call Create_Node(Tree % root % right, center+length/2., length/2.)    
  end subroutine INIT_TREE
    
end program tree_sim 
