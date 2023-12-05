program tree_sim
  implicit none
  type TreeNode
     real(8) :: Center(2)  ! Center of the line segment
     real(8) :: length           ! Length of the box
     integer :: state_vector_index !Bootleg pointer lmfao.
     real(8) :: total_mass
     real(8) :: center_mass(2)     
     type (TreeNode), pointer :: NW, NE, SW, SE
  end type TreeNode
  
  type Quadtree
     type (TreeNode), pointer :: root     
  end type Quadtree

  type (Quadtree) :: Tree
  integer, parameter :: system_size=49979
  real(8), parameter :: theta = 0.5
  real(8), parameter :: eps = 0.15 !
  !gravitational constant in (pc)^3 (Solar Mass)^-1 (seconds)^-2
  real(8), parameter :: pc_to_m = 3.2407792896664e-17 
  real(8), parameter :: Solar_to_kg = 5.0289921396853e-31
  real(8), parameter :: G =  6.674e-11 * (pc_to_m)**3/(Solar_to_kg) 
  
  real(8), parameter :: dt = 2.5e9 !Seconds
  real(8) :: SYSTEM_CENTER(2)
  !  1, 2, 3,  4,  5,  6,  7
  real(8) :: STATES(system_size, 7)  ! (X, Y, Vx, Vy, Fx, Fy, Mass)
  integer :: n, nn
  character(len=1024) :: filename
  

  SYSTEM_CENTER = (/ 0.0, 0.0 /)

  write(*,*) G
  
  !open(file='/home/alec/Backup/2d_Tree_data_/frame_ 68000_.dat', status='old', action='read', unit=12)
  open(file='initial_data.dat', status='old', action='read', unit=12)
  do n=1, system_size
     read (12, *) STATES(n, :)
     !write(*,*) STATES(n, :)
  enddo
  close(unit=12)
  write(*,*) 'System Setup'
  !call exit(-1)
  
  
  do nn = 1, 36600
     call CREATE_TREE(Tree)
     !write(*,*) 'Tree Created'
     call UPDATE_VIRTUAL_NODES(Tree)
     !write(*,*) 'Virtual Nodes updated'
     do n=1, system_size
        STATES(n, 5:6) = calculate_force(Tree, n) 
     enddo
     !write(*,*) 'Force Calculated'
     STATES(:, 1) = STATES(:, 1) + STATES(:, 3)*dt + 0.5*STATES(:, 5)/STATES(:, 7) *dt *dt
     STATES(:, 2) = STATES(:, 2) + STATES(:, 4)*dt + 0.5*STATES(:, 6)/STATES(:, 7) *dt *dt
     STATES(:, 3) = STATES(:, 3) +  dt * STATES(:, 5)/STATES(:, 7)
     STATES(:, 4) = STATES(:, 4) +  dt * STATES(:, 6)/STATES(:, 7)
     STATES(:, 5:6) = 0.0

     if (mod(nn, 100) == 0 .or. nn == 1) then

        write(filename, 42) nn
42      FORMAT('/home/alec/Backup/2d_Tree_data_/frame_', I14, '_.dat')
        
        open(file=trim(filename), unit=100)
        do n=1, system_size
           write(100, *) STATES(n, :)
        enddo
        write(*,*) nn
     endif
     
     SYSTEM_CENTER = Tree % root % center_mass
     call DESTROY_TREE(Tree)
     !write(*,*) 'Tree Destroyed'
  enddo
  close(unit=10)
contains


  real(8) function kinetic()
    integer :: i
    real(8) v
    do i = 1, system_size
       v = sqrt(STATES(i, 3)**2 + STATES(i, 4)**2)  ! Velocity magnitude
       kinetic = kinetic + 0.5 * STATES(i, 7) * v**2
    end do
  end function kinetic

  
  real(8) function potential()
    integer :: i, j
    real(8) m_i, m_j, r
    do i = 1, system_size-1
    do j = i+1, system_size
        m_i = STATES(i, 7)
        m_j = STATES(j, 7)
        r = sqrt((STATES(i, 1) - STATES(j, 1))**2 + (STATES(i, 2) - STATES(j, 2))**2)
        if (r /= 0.0) then
            potential = potential - G * m_i * m_j / r
         end if
      end do
   end do
 end function potential
  
    
  
  function grav(p1, p2, m1, m2)
    real(8), dimension(2) :: grav
    real(8) :: p1(2), p2(2), m1, m2
    real(8) :: r(2), distance
    r = p1 - p2
    distance = sqrt(r(1)*r(1)+r(2)*r(2) + eps**2)
    grav = G * m1 * m2 * r / distance**3
  end function grav
    
  
  recursive function calculate_force_(Tree, node, state_vector_index) RESULT(force)
    type (QuadTree) :: Tree
    type (TreeNode), pointer :: node
    integer :: state_vector_index
    real(8) :: dist(2), p1(2), p2(2), m1, m2
    real(8) :: dist_sqr
    real(8), dimension(2) :: force

    if (.not. associated(node)) then
       force = 0.0
    else
       
       dist = abs(node % center_mass - STATES(state_vector_index, 1:2))
       dist_sqr = sqrt(dist(1)*dist(1)+dist(2)*dist(2))

       if(state_vector_index == node % state_vector_index) then
          force = 0.0
       else
          if (.not. external_node(node) .and. node % length / dist_sqr < theta) then
             p1 = STATES(state_vector_index, 1:2)
             p2 = node % center_mass
             m1 = STATES(state_vector_index, 7)
             m2 = node % total_mass
             force = -1.0*grav(p1, p2,  m1, m2)
          else  if (external_node(node) .and. node % state_vector_index > 0) then
             p1 = STATES(state_vector_index, 1:2)
             p2 = STATES(node % state_vector_index, 1:2)
             m1 = STATES(state_vector_index, 7)
             m2 = STATES(node % state_vector_index, 7)
             force = -1.0*grav(p1, p2,  m1, m2)         
          else
             force = calculate_force_(Tree, node % NW, state_vector_index) + &
                  calculate_force_(Tree, node % NE, state_vector_index) + &
                  calculate_force_(Tree, node % SW, state_vector_index) + &
                  calculate_force_(Tree, node % SE, state_vector_index)
          endif
       endif
    endif
  end function calculate_force_
  
  recursive function calculate_force(Tree, state_vector_index) result(force)
    real(8), dimension(2) :: force
    type (QuadTree) :: Tree
    integer :: state_vector_index
    force = CALCULATE_FORCE_(Tree, Tree % root, state_vector_index)    
  end function calculate_force

  SUBROUTINE CREATE_TREE(Tree)
    type (QuadTree) :: Tree
    integer :: n
    call INIT_TREE(Tree,  (/ 0.0d0, 0.0d0 /), maxval((/ maxval(STATES(:, 1)), abs(minval(STATES(:, 1))), &
         maxval(STATES(:, 2)), abs(minval(STATES(:, 2))) /)))
    do n=1, system_size
       call INSERT(Tree, n)  
    enddo
  end SUBROUTINE CREATE_TREE
  
  RECURSIVE SUBROUTINE DESTROY_TREE_(Tree, Node)
    type (QuadTree) :: Tree
    type (TreeNode), pointer :: node
    if (associated(node)) then
       call DESTROY_TREE_(Tree, node % NW)
       call DESTROY_TREE_(Tree, node % NE)
       call DESTROY_TREE_(Tree, node % SW)
       call DESTROY_TREE_(Tree, node % SE)       
       deallocate(node)       
    endif       
  END SUBROUTINE DESTROY_TREE_
  
  RECURSIVE SUBROUTINE DESTROY_TREE(Tree)
    type (QuadTree) :: Tree
    call DESTROY_TREE_(Tree, Tree % root)
  END SUBROUTINE DESTROY_TREE
 
    
  RECURSIVE SUBROUTINE UPDATE_VIRTUAL_NODES_(Tree, node)
    type (QuadTree) :: Tree
    type (TreeNode) :: node
    real(8) :: TM, CM(2)


    if (external_node(node) .and. node % state_vector_index < 0) then
       node % total_mass = 0.0
       node % center_mass = (/ 0.0, 0.0 /)       
    else if (external_node(node) .and. node % state_vector_index > 0) then
       node % total_mass = STATES(node % state_vector_index, 7)
       node % center_mass = STATES(node % state_vector_index, 1:2)
    else if (internal_node(node)) then
       call UPDATE_VIRTUAL_NODES_(Tree, node % NW)
       call UPDATE_VIRTUAL_NODES_(Tree, node % NE)
       call UPDATE_VIRTUAL_NODES_(Tree, node % SW)
       call UPDATE_VIRTUAL_NODES_(Tree, node % SE)

       TM = node % NW % total_mass + node % NE % total_mass + &
            node % SW % total_mass + node % SE % total_mass
       CM =  (node % NW % center_mass * node % NW % total_mass + &
            node % NE % center_mass * node % NE % total_mass + &
            node % SW % center_mass * node % SW % total_mass + &
            node % SE % center_mass * node % SE % total_mass) / TM
       node % total_mass = TM
       node % center_mass = TM
    endif

  END SUBROUTINE UPDATE_VIRTUAL_NODES_
  
  RECURSIVE SUBROUTINE UPDATE_VIRTUAL_NODES(Tree)
    type (QuadTree) :: Tree
    call UPDATE_VIRTUAL_NODES_(Tree, Tree % root)
  END SUBROUTINE UPDATE_VIRTUAL_NODES

  SUBROUTINE WRITE_NODE(node)
    type (TreeNode), pointer :: node
    write(*,*) 'center', node % center, 'length', node % length, &
         'total mass', node % total_mass, 'center mass', &
         node % center_mass, 'LEAF', node % state_vector_index > 0, 'NDX', node % state_vector_index 
  END SUBROUTINE WRITE_NODE
     
  RECURSIVE SUBROUTINE WRITE_TREE(node)
    type (TreeNode), pointer :: node
    if (associated(node)) then
       call WRITE_TREE(node % NW)
       call WRITE_TREE(node % NE)       
       call WRITE_NODE(node)
       call WRITE_TREE(node % SW)
       call WRITE_TREE(node % SE)      
    endif    
  end SUBROUTINE WRITE_TREE
  
  RECURSIVE SUBROUTINE WRITE_LEAVES(node)
    type (TreeNode), pointer :: node
    if (associated(node)) then
       call WRITE_TREE(node % NW)
       call WRITE_TREE(node % NE)       
       if (node % state_vector_index > 0) then
          write(*,*) STATES(node % state_vector_index, :)
       endif
       call WRITE_TREE(node % SW)
       call WRITE_TREE(node % SE)
    endif
  end SUBROUTINE WRITE_LEAVES


  logical function internal_node(Node)
    type (TreeNode) :: Node
    internal_node = associated(Node % NW) .or. associated(Node % NE) .or. &
         associated(Node % SW) .or. associated(Node % SE)
  end function internal_node

  
  logical function external_node(Node)
    type (TreeNode) :: Node
    external_node = .not. associated(Node % NW) .and. .not. associated(Node % NE) .and. &
         .not. associated(Node % SW) .and. .not. associated(Node % SE)
  end function external_node

  RECURSIVE subroutine INTERNAL_INSERT_(Tree, current, state_vector_index)
    type (QuadTree) :: Tree
    type (TreeNode) :: current
    integer :: state_vector_index

    if (current % center(1) <= STATES(state_vector_index, 1) .and. &
         current % center(2) <= STATES(state_vector_index, 2)) then
       call INSERT_(Tree, current % NE, state_vector_index)
    else if (current % center(1) <= STATES(state_vector_index, 1) .and. &
         current % center(2) >= STATES(state_vector_index, 2)) then
       call INSERT_(Tree, current % SE, state_vector_index)
    else if (current % center(1) >= STATES(state_vector_index, 1) .and. &
         current % center(2) >= STATES(state_vector_index, 2)) then
       call INSERT_(Tree, current % SW, state_vector_index)
    else if (current % center(1) >= STATES(state_vector_index, 1) .and. &
         current % center(2) <= STATES(state_vector_index, 2)) then
       call INSERT_(Tree, current % NW, state_vector_index)
    !else
       !write(*,*) 'Uncaught Case :( '
       !write(*,*) state_vector_index, STATES(state_vector_index, 1:2)
    endif
       
       
  end subroutine INTERNAL_INSERT_

  RECURSIVE subroutine EXTERNAL_INSERT_(Tree, current, state_vector_index)
    type (QuadTree) :: Tree
    type (TreeNode) :: current
    integer :: state_vector_index
    integer :: old
    real(8) :: center(2), length
    
    old = current % state_vector_index
    current % state_vector_index = -1 !(/ 0.0, 0.0, 0.0, -1.0 /)
    center = current % center
    length = current % length
    
    call Create_Node(current % NE, center+length/2., length/2.)
    call Create_Node(current % SW, center-length/2., length/2.)
    call Create_Node(current % NW, (/ center(1) -length/2., center(2) + length/2 /), length/2.)
    call Create_Node(current % SE, (/ center(1) +length/2., center(2) - length/2 /), length/2.)

    
    call INSERT_(Tree, Tree % root, state_vector_index)
    call INSERT_(Tree, Tree % root, old)    
  end subroutine EXTERNAL_INSERT_

  RECURSIVE subroutine INSERT_(Tree, current, state_vector_index)
    type (QuadTree) :: Tree
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
    !else
     !  write(*,*) 'Uncaught case in insert'
    endif
  end subroutine INSERT_

  subroutine INSERT(Tree, state_vector_index)
    type (QuadTree) :: Tree
    integer :: state_vector_index
    !real(8) :: top_bound_left, top_bound_right, bot_bound_left, bot_bound_right

    !top_bound_left = Tree % root % center(1) - Tree % root % length
    !top_bound_right = Tree % root % center(1) + Tree % root % length

    !bot_bound_left = Tree % root % center(2) - Tree % root % length
    !bot_bound_right = Tree % root % center(2) + Tree % root % length

    
    
    !if (STATES(state_vector_index, 1) < top_bound_left .or. &
    !STATES(state_vector_index, 1) > top_bound_right .or. &
    !STATES(state_vector_index, 2) < bot_bound_left .or. &
    !STATES(state_vector_index, 2) > bot_bound_right ) then
    !write(*,*) 'point outside of tree bounds... Exiting.'
    !call exit(-1)
    !else
    call INSERT_(Tree, Tree % root, state_vector_index)
    !endif
    
  end subroutine INSERT
  
  subroutine CREATE_NODE(node, center, length)
    type (TreeNode), pointer :: node
    real(8) :: center(2), length
    allocate (node)
    node % center = center
    node % length = length
    node % state_vector_index = -1
    nullify (node % NW)
    nullify (node % NE)
    nullify (node % SW)
    nullify (node % SE)
  end subroutine CREATE_NODE
  
  subroutine INIT_TREE(Tree, center, length)
    real(8) :: center(2), length
    type (QuadTree) :: Tree
    call Create_Node(Tree % root, center, length)

    call Create_Node(Tree % root % NE, center+length/2., length/2.)
    call Create_Node(Tree % root % SW, center-length/2., length/2.)
    call Create_Node(Tree % root % NW, (/ center(1) -length/2., center(2) + length/2 /), length/2.)
    call Create_Node(Tree % root % SE, (/ center(1) +length/2., center(2) - length/2 /), length/2.)
  
  end subroutine INIT_TREE
    
end program tree_sim 
