!dec$ if .not. defined(__ParticleSimu_general_sub_f90)
!dec$ define __ParticleSimu_general_f90

module ParticleSimu_general_sub
    implicit none
    contains
    
        !****************************************************************************
        !  sort the array according to the given order
        !  array        :   real(kind=8),allocatable, array to be sorted
        !  dimension    :   integer(kind=4), indicate which column of array should be sorted
        !  sequence     :   character(len=1), 'u' for up order, 'd' for down order
        !****************************************************************************
        subroutine sort(array,dimension,sequence)
            implicit none
            real(kind=8),allocatable    ::  array(:,:)
            integer(kind=4),intent(in)  ::  dimension
            character(len=1),intent(in) ::  sequence
            integer(kind=4)             ::  i,j,k
            integer(kind=4)             ::  row_size, column_size
            real(kind=8)                ::  temp
            
            row_size = size(array,1)
            column_size = size(array,2)
            if (sequence .eq. 'u') then
                do i = 1, row_size - 1
                    k = i
                    do j = i + 1, row_size
                        if (array(k,dimension) .gt. array(j,dimension)) then
                            k = j
                        end if
                    end do
                    if (k .ne. i) then
                        do j = 1, column_size
                            temp = array(k,j)
                            array(k,j) = array(i,j)
                            array(i,j) = array(k,j)
                        end do
                    end if
                end do
            end if

            if (sequence .eq. 'd') then
                do i = 1, row_size - 1
                    k = i
                    do j = i + 1, row_size
                        if (array(k,dimension) .lt. array(j,dimension)) then
                            k = j
                        end if
                    end do
                    if (k .ne. i) then
                        do j = 1, column_size
                            temp = array(k,j)
                            array(k,j) = array(i,j)
                            array(i,j) = array(k,j)
                        end do
                    end if
                end do
            end if
        end subroutine

        !****************************************************************************
        !
        !  calculate the area of polygon
        !  polygon_vertex   size:   [2,total_vertex_num]
        !
        !****************************************************************************
        function polygon_area(polygon_vertex)
            implicit none
            real(kind=8),allocatable    ::  polygon_vertex(:,:)
            integer(kind=4)             ::  total_vertex_num, i, ii, jj
            real(kind=8)                ::  polygon_area, sum, temp(2,3)

            total_vertex_num = size(polygon_vertex,2)
            temp(1,1) = polygon_vertex(1,1)
            temp(2,1) = polygon_vertex(2,1)
            sum = 0.0
            do i = 2, total_vertex_num-1
                temp(1,2) = polygon_vertex(1,i)
                temp(2,2) = polygon_vertex(2,i)
                temp(1,3) = polygon_vertex(1,i+1)
                temp(2,3) = polygon_vertex(2,i+1)
                sum = sum + triangle_area(temp)
            end do

            polygon_area = sum

            return

        end function polygon_area

        !****************************************************************************
        !
        !  calculate the area of triangle
        !  triangle_vertex  size:   [2,3]
        !
        !****************************************************************************
        function triangle_area(triangle_vertex)
            implicit none
            real(kind=8)    ::  triangle_vertex(2,3)
            real(kind=8)    ::  triangle_area
            real(kind=8)    ::  x1, x2, x3, y1, y2, y3

            x1 = triangle_vertex(1,1)
            x2 = triangle_vertex(1,2)
            x3 = triangle_vertex(1,3)
            y1 = triangle_vertex(2,1)
            y2 = triangle_vertex(2,2)
            y3 = triangle_vertex(2,3)
            triangle_area = 0.5 * (-x2 * y1 + x3 * y1 + x1 * y2 - x3 * y2 - x1 * y3 + x2 * y3)

            return

        end function triangle_area

        !****************************************************************************
        !
        !  determine whether two line segments are crossed (on 2d plane)
        !  x (y)    size:   [4] x1, x2, x3, x4
        !
        !****************************************************************************
        function is_line_cross(x,y)
            implicit none
            real(kind=8)    ::  x(4), y(4)
            logical         ::  is_line_cross
            real(kind=8)    ::  delta, s, t
            
            delta = (x(3) - x(4)) * (y(2) - y(1)) - (x(2) - x(1)) * (y(3) - y(4))
            if (delta .ne. 0.0) then
                s = ((x(3) - x(1)) * (y(4) - y(3)) - (x(4) - x(3)) * (y(3) - y(1))) / delta
                t = ((x(3) - x(1)) * (y(2) - y(1)) - (x(2) - x(1)) * (y(3) - y(1))) / delta
                if (((s .ge. 0.0) .and. (s .le. 1.0)) .and. ((t .ge. 0.0) .and. (t .le. 1.0))) then
                    is_line_cross = .true.
                    return
                end if
            end if
            is_line_cross = .false.
            return
        end function is_line_cross
        
        ! return the aspect ratio of polygon
        ! polygon_vertex    size:   [2,total_vertex_num]
        function aspect_ratio_of_polygon(polygon_vertex) 
            implicit none
            real(kind=8),allocatable    ::  polygon_vertex(:,:)
            real(kind=8),allocatable    ::  aspect_ratio_for_each_side(:), perpendicular_distance(:), foot_point(:,:)
            real(kind=8)                ::  x1, x2, x3, y1, y2, y3
            real(kind=8)                ::  cross_x_pos, cross_y_pos, max_perp_dis, max_hori_dis, hori_dis
            real(kind=8)                ::  aspect_ratio_of_polygon
            integer(kind=4)             ::  i, j, ii, jj, side_num
            
            side_num = size(polygon_vertex,2)
            if (allocated(aspect_ratio_for_each_side)) deallocate(aspect_ratio_for_each_side)
            allocate(aspect_ratio_for_each_side(side_num))
            if (allocated(perpendicular_distance)) deallocate(perpendicular_distance)
            allocate(perpendicular_distance(side_num))
            if (allocated(foot_point)) deallocate(foot_point)
            allocate(foot_point(2,side_num))
            
            do i = 1, side_num-1
                x1 = polygon_vertex(1,i)
                y1 = polygon_vertex(2,i)
                x2 = polygon_vertex(1,i+1)
                y2 = polygon_vertex(2,i+1)
                do j = 1, side_num
                    if ((j .ne. i) .and. (j .ne. i+1)) then
                        x3 = polygon_vertex(1,j)
                        y3 = polygon_vertex(2,j)
                        cross_x_pos = ((x2 - x1) ** 2.0 * x3 + (y2 - y1) ** 2.0 * x2 + (y3 - y2) * (x2 - x1) * (y2 - y1)) / ((y2 - y1) ** 2.0 + (x2 - x1) ** 2.0)
                        if (y2 .eq. y1) cross_y_pos = y1
                        if (y2 .ne. y1) cross_y_pos = -(x2 - x1) / (y2 - y1) * (cross_x_pos - x3) + y3
                        perpendicular_distance(j) = ((cross_y_pos - y3) ** 2.0 + (cross_x_pos - x3) ** 2.0) ** 0.5
                        foot_point(1,j) = cross_x_pos
                        foot_point(2,j) = cross_y_pos
                    end if
                    if (j .eq. i) then
                        perpendicular_distance(j) = 0
                        foot_point(1,j) = x1
                        foot_point(2,j) = y1
                    end if
                    if (j .eq. i+1) then
                        perpendicular_distance(j)=0
                        foot_point(1,j)=x2
                        foot_point(2,j)=y2
                    end if
                end do
                max_perp_dis = maxval(perpendicular_distance)
                max_hori_dis = 0
                do ii = 1, side_num-1
                    x1 = foot_point(1,ii)
                    y1 = foot_point(2,ii)
                    do jj = ii+1, side_num
                        x2 = foot_point(1,jj)
                        y2 = foot_point(2,jj)
                        hori_dis = ((y2 - y1) ** 2.0 + (x2 - x1) ** 2.0) ** 0.5
                        if (max_hori_dis .lt. hori_dis) max_hori_dis = hori_dis
                    enddo
                enddo
                aspect_ratio_for_each_side(i) = max_perp_dis / max_hori_dis
            enddo
            
            ! special case
            x1 = polygon_vertex(1,side_num)
            y1 = polygon_vertex(2,side_num)
            x2 = polygon_vertex(1,1)
            y2 = polygon_vertex(2,1)
            perpendicular_distance = 0.0
            foot_point = 0.0
            do j = 1, side_num
                if ((j .ne. side_num) .and. (j .ne. 1)) then
                    x3 = polygon_vertex(1,j)
                    y3 = polygon_vertex(2,j)
                    cross_x_pos = ((x2 - x1) ** 2.0 * x3 + (y2 - y1) ** 2.0 * x2 + (y3 - y2) * (x2-x1) * (y2-y1)) / ((y2-y1) ** 2.0 + (x2 - x1) ** 2.0)
                    if (y2 .eq. y1) cross_y_pos = y1
                    if (y2 .ne. y1) cross_y_pos = -(x2 - x1) / (y2 - y1) * (cross_x_pos - x3) + y3
                    perpendicular_distance(j) = ((cross_y_pos - y3) ** 2.0 + (cross_x_pos - x3) ** 2.0) ** 0.5
                    foot_point(1,j) = cross_x_pos
                    foot_point(2,j) = cross_y_pos
                end if
                if (j .eq. side_num) then
                    perpendicular_distance(j) = 0
                    foot_point(1,j) = x1
                    foot_point(2,j) = y1
                end if
                if (j .eq. 1) then
                    perpendicular_distance(j) = 0
                    foot_point(1,j) = x2
                    foot_point(2,j) = y2
                end if
            end do
            max_perp_dis = maxval(perpendicular_distance)
            
            max_hori_dis = 0
            do ii = 1, side_num-1
                x1 = foot_point(1,ii);
                y1 = foot_point(2,ii);
                do jj = ii+1, side_num
                    x2 = foot_point(1,jj);
                    y2 = foot_point(2,jj);
                    hori_dis = ((y2 - y1) ** 2.0 + (x2 - x1) ** 2.0) ** 0.5
                    if (max_hori_dis .lt. hori_dis) max_hori_dis = hori_dis
                end do
            end do
            aspect_ratio_for_each_side(side_num) = max_perp_dis / max_hori_dis
            aspect_ratio_of_polygon = maxval(aspect_ratio_for_each_side)
            return
        end function

        ! extend the polygon according to the aspect_ratio (length-width ratio)
        ! CAUTION: will change the value of polygon_vertex
        ! polygon_vertex    size:   [2,total_vertex_num]
        subroutine extend_polygon(polygon_vertex,required_aspect_ratio)
            implicit none
            real(kind=8),allocatable    ::  polygon_vertex(:,:)
            real(kind=8),intent(in)     ::  required_aspect_ratio
            real(kind=8),allocatable    ::  polygon_vertex_in_ref(:,:)
            real(kind=8),allocatable    ::  aspect_ratio_for_each_side(:), perpendicular_distance(:), foot_point(:,:)
            real(kind=8)                ::  x1, x2, x3, y1, y2, y3
            real(kind=8)                ::  cross_x_pos, cross_y_pos, max_perp_dis, max_hori_dis, hori_dis
            real(kind=8)                ::  max_aspect_ratio, angle
            real(kind=8)                ::  zoom_matrix(3,3), rotation_matrix(3,3), point_in_ref(3)
            integer(kind=4)             ::  i, j, ii, jj, side_num, index
            
            side_num=size(polygon_vertex,2);
            if (allocated(polygon_vertex_in_ref)) deallocate(polygon_vertex_in_ref)
            allocate(polygon_vertex_in_ref(2,side_num))
            if (allocated(aspect_ratio_for_each_side)) deallocate(aspect_ratio_for_each_side)
            allocate(aspect_ratio_for_each_side(side_num))
            if (allocated(perpendicular_distance)) deallocate(perpendicular_distance)
            allocate(perpendicular_distance(side_num))
            if (allocated(foot_point)) deallocate(foot_point)
            allocate(foot_point(2,side_num))
            
            do i = 1, side_num-1
                x1 = polygon_vertex(1,i)
                y1 = polygon_vertex(2,i)
                x2 = polygon_vertex(1,i+1)
                y2 = polygon_vertex(2,i+1)
                do j = 1, side_num
                    if ((j .ne. i) .and. (j .ne. i+1)) then
                        x3 = polygon_vertex(1,j)
                        y3 = polygon_vertex(2,j)
                        cross_x_pos = ((x2 - x1) ** 2.0 * x3 + (y2 - y1) ** 2.0 * x2 + (y3 - y2) * (x2 - x1) * (y2 - y1)) / ((y2 - y1) ** 2.0 + (x2 - x1) ** 2.0)
                        if (y2 .eq. y1) cross_y_pos = y1
                        if (y2 .ne. y1) cross_y_pos = -(x2 - x1) / (y2 - y1) * (cross_x_pos - x3) + y3
                        perpendicular_distance(j) = ((cross_y_pos - y3) ** 2.0 + (cross_x_pos - x3) ** 2.0) ** 0.5
                        foot_point(1,j) = cross_x_pos
                        foot_point(2,j) = cross_y_pos
                    end if
                    if (j .eq. i) then
                        perpendicular_distance(j) = 0
                        foot_point(1,j) = x1
                        foot_point(2,j) = y1
                    end if
                    if (j .eq. i+1) then
                        perpendicular_distance(j) = 0
                        foot_point(1,j) = x2
                        foot_point(2,j) = y2
                    end if
                end do
                max_perp_dis = maxval(perpendicular_distance)
                max_hori_dis = 0
                do ii = 1, side_num-1
                    x1 = foot_point(1,ii)
                    y1 = foot_point(2,ii)
                    do jj = ii+1, side_num
                        x2 = foot_point(1,jj)
                        y2 = foot_point(2,jj)
                        hori_dis = ((y2 - y1) ** 2.0 + (x2 - x1) ** 2.0) ** 0.5
                        if (max_hori_dis .lt. hori_dis) max_hori_dis = hori_dis
                    end do
                end do
                aspect_ratio_for_each_side(i) = max_perp_dis / max_hori_dis
            end do
            
            ! the case of i=side_num (last point and first point)
            x1 = polygon_vertex(1,side_num)
            y1 = polygon_vertex(2,side_num)
            x2 = polygon_vertex(1,1)
            y2 = polygon_vertex(2,1)
            perpendicular_distance = 0.0
            foot_point = 0.0
            do j = 1, side_num
                if ((j .ne. side_num) .and. (j .ne. 1)) then
                    x3 = polygon_vertex(1,j)
                    y3 = polygon_vertex(2,j)
                    cross_x_pos = ((x2 - x1) ** 2.0 * x3 + (y2 - y1) ** 2.0 * x2 + (y3 - y2) * (x2 - x1) * (y2 - y1)) / ((y2 - y1) ** 2.0 + (x2 - x1) ** 2.0)
                    if (y2 .eq. y1) cross_y_pos = y1
                    if (y2 .ne. y1) cross_y_pos = -(x2 - x1) / (y2 - y1) * (cross_x_pos - x3) + y3
                    perpendicular_distance(j) = ((cross_y_pos - y3) ** 2.0 + (cross_x_pos - x3) ** 2.0) ** 0.5
                    foot_point(1,j) = cross_x_pos
                    foot_point(2,j) = cross_y_pos
                end if
                if (j .eq. side_num) then
                    perpendicular_distance(j) = 0
                    foot_point(1,j) = x1
                    foot_point(2,j) = y1
                end if
                if (j .eq. 1) then
                    perpendicular_distance(j) = 0
                    foot_point(1,j) = x2
                    foot_point(2,j) = y2
                end if
            end do
            max_perp_dis = maxval(perpendicular_distance)
            
            max_hori_dis = 0
            do ii = 1, side_num-1
                x1 = foot_point(1,ii)
                y1 = foot_point(2,ii)
                do jj = ii+1, side_num
                    x2 = foot_point(1,jj)
                    y2 = foot_point(2,jj)
                    hori_dis = ((y2 - y1) ** 2.0 + (x2 - x1) ** 2.0) ** 0.5
                    if (max_hori_dis .lt. hori_dis) max_hori_dis = hori_dis
                end do
            end do
            aspect_ratio_for_each_side(side_num) = max_perp_dis / max_hori_dis
            
            ! find the max. aspect ratio and its corresponding side
            index = 1
            max_aspect_ratio = aspect_ratio_for_each_side(1)
            do i = 2, side_num
                if (max_aspect_ratio .lt. aspect_ratio_for_each_side(i)) then
                    index = i
                    max_aspect_ratio = aspect_ratio_for_each_side(i)
                end if
            end do
            
            polygon_vertex_in_ref = polygon_vertex
            
            ! Find the points of ref. side according to index
            if (index .lt. side_num) then
                x1 = polygon_vertex(1,index)      ! First point of ref. side
                y1 = polygon_vertex(2,index)
                x2 = polygon_vertex(1,index+1)    ! Second point of ref. side
                y2 = polygon_vertex(2,index+1)
            else
                x1 = polygon_vertex(1,side_num)
                y1 = polygon_vertex(2,side_num)
                x2 = polygon_vertex(1,1)
                y2 = polygon_vertex(2,1)
            endif
            
            angle = atan((y2 - y1) / (x2 - x1))
            
            ! global pos to local pos
            zoom_matrix = 0.0
            zoom_matrix(1,1) = 1.0
            zoom_matrix(2,2) = 1.0
            zoom_matrix(3,3) = 1.0
            zoom_matrix(3,1) = -x1
            zoom_matrix(3,2) = -y1
            
            rotation_matrix = 0.0
            rotation_matrix(1,1) = cos(-angle)
            rotation_matrix(1,2) = sin(-angle)
            rotation_matrix(2,1) =-sin(-angle)
            rotation_matrix(2,2) = cos(-angle)
            rotation_matrix(3,3) = 1.0
            
            do i = 1, side_num
                point_in_ref = 0.0
                point_in_ref(1) = polygon_vertex(1,i)
                point_in_ref(2) = polygon_vertex(2,i)
                point_in_ref(3) = 1.0
                point_in_ref = matmul(matmul(point_in_ref,zoom_matrix),rotation_matrix)
                polygon_vertex_in_ref(1,i) = point_in_ref(1)
                polygon_vertex_in_ref(2,i) = point_in_ref(2)
            end do
            do i = 1, side_num
                polygon_vertex_in_ref(2,i) = polygon_vertex_in_ref(2,i) * required_aspect_ratio / max_aspect_ratio
            end do

            ! back to global coordinate
            zoom_matrix = 0.0
            zoom_matrix(1,1) = 1.0
            zoom_matrix(2,2) = 1.0
            zoom_matrix(3,3) = 1.0
            zoom_matrix(3,1) = x1
            zoom_matrix(3,2) = y1
            
            rotation_matrix = 0.0
            rotation_matrix(1,1) = cos(angle)
            rotation_matrix(1,2) = sin(angle)
            rotation_matrix(2,1) =-sin(angle)
            rotation_matrix(2,2) = cos(angle)
            rotation_matrix(3,3) = 1.0
            
            do i = 1, side_num
                point_in_ref = 0.0
                point_in_ref(1) = polygon_vertex_in_ref(1,i)
                point_in_ref(2) = polygon_vertex_in_ref(2,i)
                point_in_ref(3) = 1.0
                point_in_ref = matmul(matmul(point_in_ref,zoom_matrix),rotation_matrix)
                polygon_vertex_in_ref(1,i) = point_in_ref(1)
                polygon_vertex_in_ref(2,i) = point_in_ref(2)
            end do
            polygon_vertex = polygon_vertex_in_ref
        end subroutine

        ! expand the polygon according to the given area
        ! the shape of polygon remains the same
        ! CAUTION: will change the value of polygon_vertex
        ! polygon_vertex    size:   [2,total_vertex_num]
        subroutine expand_polygon(polygon_vertex,area)
            implicit none
            real(kind=8),allocatable    ::  polygon_vertex(:,:)
            real(kind=8),intent(in)     ::  area
            real(kind=8),allocatable    ::  expanded_polygon_vertex(:,:)
            integer(kind=4)             ::  i, side_num
            real(kind=8)                ::  origin_area, linear_scale_ratio, k, distance, x_temp_pos1, x_temp_pos2
            
            side_num = size(polygon_vertex,2)
            origin_area = polygon_area(polygon_vertex)
            linear_scale_ratio = (area / origin_area) ** 0.5
            
            if (allocated(expanded_polygon_vertex)) deallocate(expanded_polygon_vertex)
            allocate(expanded_polygon_vertex(2,side_num))
            expanded_polygon_vertex(1,1) = polygon_vertex(1,1)
            expanded_polygon_vertex(2,1) = polygon_vertex(2,1)
            
            do i = 2, side_num
                k = (polygon_vertex(2,i) - expanded_polygon_vertex(2,1)) / (polygon_vertex(1,i) - expanded_polygon_vertex(1,1))
                distance = linear_scale_ratio * (((polygon_vertex(2,i) - expanded_polygon_vertex(2,1)) ** 2.0 + (polygon_vertex(1,i) - expanded_polygon_vertex(1,1)) ** 2.0) ** 0.5)
                x_temp_pos1 =  distance / ((1.0 + k ** 2.0) ** 0.5) + expanded_polygon_vertex(1,1)
                x_temp_pos2 = -distance / ((1.0 + k ** 2.0) ** 0.5) + expanded_polygon_vertex(1,1)
                if (linear_scale_ratio .ge. 1.0) then
                    if ((polygon_vertex(1,i) - x_temp_pos1) * (polygon_vertex(1,i) - expanded_polygon_vertex(1,1)) .lt. 0.0) then
                        expanded_polygon_vertex(1,i) = x_temp_pos1
                    else
                        expanded_polygon_vertex(1,i) = x_temp_pos2
                    end if
                    expanded_polygon_vertex(2,i) = k * (expanded_polygon_vertex(1,i) - expanded_polygon_vertex(1,1)) + expanded_polygon_vertex(2,1)
                end if
                if (linear_scale_ratio .lt. 1.0) then
                    if ((x_temp_pos1 - polygon_vertex(1,i)) * (x_temp_pos1 - expanded_polygon_vertex(1,1)) .lt. 0.0) then
                        expanded_polygon_vertex(1,i) = x_temp_pos1
                    else
                        expanded_polygon_vertex(1,i) = x_temp_pos2
                    end if
                    expanded_polygon_vertex(2,i) = k * (expanded_polygon_vertex(1,i) - expanded_polygon_vertex(1,1)) + expanded_polygon_vertex(2,1)
                end if
            end do
            polygon_vertex = expanded_polygon_vertex
        end subroutine expand_polygon

        !****************************************************************************
        !
        !  determine whether a point in a polygon
        !  polygon_vertex   size:   [2,total_vertex_num]
        !  point            size:   [2]
        !
        !****************************************************************************
        function is_point_in_polygon(polygon_vertex,point)
            implicit none
            real(kind=8),allocatable    ::  polygon_vertex(:,:)
            real(kind=8)                ::  point(2)
            integer(kind=4)             ::  total_side_num, i
            real(kind=8)                ::  delta, s, t
            real(kind=8)                ::  polygon_origin_area, polygon_new_area
            real(kind=8)                ::  triangle_pos(2,3)
            logical                     ::  is_point_in_polygon
            
            total_side_num = size(polygon_vertex,2)
            
            polygon_origin_area = polygon_area(polygon_vertex)
            polygon_new_area = 0.0
            
            triangle_pos(1,1) = point(1)
            triangle_pos(2,1) = point(2)
            do i = 1, total_side_num-1
                triangle_pos(1,2) = polygon_vertex(1,i)
                triangle_pos(2,2) = polygon_vertex(2,i)
                triangle_pos(1,3) = polygon_vertex(1,i+1)
                triangle_pos(2,3) = polygon_vertex(2,i+1)
                polygon_new_area = polygon_new_area + abs(triangle_area(triangle_pos))
            enddo
            triangle_pos(1,2) = polygon_vertex(1,total_side_num)
            triangle_pos(2,2) = polygon_vertex(2,total_side_num)
            triangle_pos(1,3) = polygon_vertex(1,1)
            triangle_pos(2,3) = polygon_vertex(2,1)
            polygon_new_area = polygon_new_area + abs(triangle_area(triangle_pos))
            
            if (abs((polygon_new_area-polygon_origin_area)/polygon_origin_area) .lt. 0.001) then            
                is_point_in_polygon = .true.
            else
                is_point_in_polygon = .false.
            end if
            return

        end function is_point_in_polygon

        !******************************************************************************************
        !
        !  determine whether two polygons are overlapped with each other 
        !  polygon_a_vertex(polygon_b_vertex)       :   [2,total_a_vertex_num(total_b_vertex_num)]
        !
        !******************************************************************************************
        function is_polygon_overlapped(polygon_a_vertex,polygon_b_vertex)
            implicit none
            real(kind=8),allocatable    ::  polygon_a_vertex(:,:), polygon_b_vertex(:,:)
            integer(kind=4)             ::  i,j
            real(kind=8)                ::  vertex_x_pos,vertex_y_pos
            integer(kind=4)             ::  is_polygon_overlapping
            real(kind=8)                ::  point(2), x(4), y(4)
            integer(kind=4)             ::  polygon_a_side_num, polygon_b_side_num
            logical                     ::  is_polygon_overlapped
            
            polygon_a_side_num = size(polygon_a_vertex,2)
            polygon_b_side_num = size(polygon_b_vertex,2)

            ! determine whether a is in b
            do i = 1, polygon_a_side_num
                point(1) = polygon_a_vertex(1,i)
                point(2) = polygon_a_vertex(2,i)
                if (is_point_in_polygon(polygon_b_vertex,point) .eq. .true.) then
                    is_polygon_overlapped = .true.
                    return
                end if
            end do
            
            ! determine whether b is in a
            do i = 1, polygon_b_side_num
                point(1) = polygon_b_vertex(1,i)
                point(2) = polygon_b_vertex(2,i)
                if (is_point_in_polygon(polygon_a_vertex,point) .eq. .true.) then
                    is_polygon_overlapped = .true.
                    return
                end if
            end do
            
            ! detect whether any lines are intersected
            x = 0.0
            y = 0.0
            do i = 1, polygon_a_side_num-1
                x(1) = polygon_a_vertex(1,i)
                x(2) = polygon_a_vertex(1,i+1)
                y(1) = polygon_a_vertex(2,i)
                y(2) = polygon_a_vertex(2,i+1)
                do j = 1, polygon_b_side_num-1
                    x(3) = polygon_b_vertex(1,j)
                    x(4) = polygon_b_vertex(1,j+1)
                    y(3) = polygon_b_vertex(2,j)
                    y(4) = polygon_b_vertex(2,j+1)
                    if (is_line_cross(x,y) .eq. .true.) then
                        is_polygon_overlapped = .true.
                        return
                    end if
                end do
                x(3) = polygon_b_vertex(1,polygon_b_side_num)
                x(4) = polygon_b_vertex(1,1)
                y(3) = polygon_b_vertex(2,polygon_b_side_num)
                y(4) = polygon_b_vertex(2,1)
                if (is_line_cross(x,y) .eq. .true.) then
                    is_polygon_overlapped = .true.
                    return
                end if
            end do
            
            x(1) = polygon_a_vertex(1,polygon_a_side_num)
            x(2) = polygon_a_vertex(1,1)
            y(1) = polygon_a_vertex(2,polygon_a_side_num)
            y(2) = polygon_a_vertex(2,1)
            do j = 1, polygon_b_side_num-1
                x(3) = polygon_b_vertex(1,j)
                x(4) = polygon_b_vertex(1,j+1)
                y(3) = polygon_b_vertex(2,j)
                y(4) = polygon_b_vertex(2,j+1)
                if (is_line_cross(x,y) .eq. .true.) then
                    is_polygon_overlapped = .true.
                    return
                endif
            enddo
            x(3) = polygon_b_vertex(1,polygon_b_side_num)
            x(4) = polygon_b_vertex(1,1)
            y(3) = polygon_b_vertex(2,polygon_b_side_num)
            y(4) = polygon_b_vertex(2,1)
            if (is_line_cross(x,y) .eq. .true.) then
                is_polygon_overlapped = .true.
                return
            endif
            is_polygon_overlapped = .false.
            return
        end function is_polygon_overlapped

        !******************************************************************************************
        !
        !  determine whether two ellipses are overlapped with each other 
        !  ellipse_a (ellipse_b)       :   [5]      radius_a, radius_b, angle, center_x_pos, center_y_pos
        !
        !  use the CHARACTERISTIC POLYNOMIAL to conduct the seperation check
        !  in analytical geometry, an arbitrary ellipse is defined as transpose(X)AX=0 where A is
        !  the symmetric matrix of coefficients
        !  Assuming A is denoted as
        !  | A      B/2     D/2 |
        !  | B/2    C       E/2 |
        !  | D/2    E/2     F |
        !  The equation of the ellipse will be :
        !  A X^2 + B X Y + C Y^2 + D X + E Y + F = 0
        !******************************************************************************************
        function is_ellipse_overlapped(ellipse_a,ellipse_b)
            implicit none
            real(kind=8),intent(in)     ::  ellipse_a(5), ellipse_b(5)
            real(kind=8)                ::  coefficient(6)
            real(kind=8)                ::  A1, B1, C1, D1, E1, F1, A2, B2, C2, D2, E2, F2
            real(kind=8)                ::  A, B, C, para_A, para_B, para_C, para_D, coeff1, coeff2, coeff3
            logical                     ::  is_ellipse_overlapped

            call get_ellipse_general_equation_coefficient(ellipse_a(1),ellipse_a(2),ellipse_a(3),ellipse_a(4),ellipse_a(5),coefficient)

            ! A*X^2 + B*X*Y + C*Y^2 + D*X + E*Y + F = 0
            A1 = coefficient(1)
            B1 = coefficient(2)
            C1 = coefficient(3)
            D1 = coefficient(4)
            E1 = coefficient(5)
            F1 = coefficient(6)

            call get_ellipse_general_equation_coefficient(ellipse_b(1),ellipse_b(2),ellipse_b(3),ellipse_b(4),ellipse_b(5),coefficient)

            A2 = coefficient(1)
            B2 = coefficient(2)
            C2 = coefficient(3)
            D2 = coefficient(4)
            E2 = coefficient(5)
            F2 = coefficient(6)

            ! form the initial characteristic polynomial
            ! f(x) = A*x^3 + B*x^2 + C*x + D
            para_A = -A1*E1**2.0/4.0 - C1*D1**2.0/4.0 - B1**2.0*F1/4.0 + A1*C1*F1 + B1*D1*E1/4.0
            para_B = -A2*E1**2.0/4.0 - C2*D1**2.0/4.0 - B1**2.0*F2/4.0 + A1*C1*F2 + A1*C2*F1 + &
                    & A2*C1*F1 - B1*B2*F1/2.0 - A1*E1*E2/2.0 + B1*D1*E2/4.0 + B1*D2*E1/4.0 + B2*D1*E1/4.0 - C1*D1*D2/2.0
            para_C = -A1*E2**2.0/4.0 - C1*D2**2.0/4.0 - B2**2.0*F1/4.0 + A1*C2*F2 + A2*C1*F2 + A2*C2*F1 - B1*B2*F2/2.0 - &
                    & A2*E1*E2/2.0 + B1*D2*E2/4.0 + B2*D1*E2/4.0 + B2*D2*E1/4.0 - C2*D1*D2/2.0 
            para_D = A2*C2*F2 - C2*D2**2.0/4.0 - B2**2.0*F2/4.0 - A2*E2**2.0/4.0  + B2*D2*E2/4.0

            ! turned monic
            ! f(x) = x^3 + a*x^2 + b*x +c
            a = para_B / para_A
            b = para_C / para_A
            c = para_D / para_A

            coeff1 = -3.0 * b + a ** 2.0
            coeff2 = 3.0 * a * c + b * a ** 2.0 - 4.0 * b ** 2.0
            coeff3 = -27 * c ** 2.0 + 18 * c * a * b + a ** 2.0 * b ** 2.0 - 4 * a ** 3.0 * c - 4 * b ** 3.0

            if ((a .ge. 0) .and. (coeff1 .gt. 0) .and. (coeff2 .lt. 0) .and. (coeff3 .gt. 0)) then
                is_ellipse_overlapped = .false.
                return
            end if

            if ((a .lt. 0) .and. (coeff1 .gt. 0) .and. (coeff3 .gt. 0)) then
                is_ellipse_overlapped = .false.
                return
            end if

            is_ellipse_overlapped = .true.
            return

        end function

        ! this subroutine calc. the coefficients of the general equation for an arbitrary ellipse
        ! the general euqation of an ellipse is like
        ! A * X^2 + B * X * Y + C * Y^2 + D * X + E * Y + F = 0
        ! a, b, angle, xc, yc   :   semi-major radius, semi-minor radius, rotation angle, center coordinates
        ! coefficient           :   [6]     A, B, C, D, E, F
        subroutine get_ellipse_general_equation_coefficient(a,b,angle,xc,yc,coefficient)
            implicit none
            real(kind=8),intent(in)     ::  a, b, angle, xc, yc
            real(kind=8),intent(out)    ::  coefficient(6)

            coefficient(1) = (a ** 2.0) * (sin(angle) ** 2.0) + (b ** 2.0) * (cos(angle) ** 2.0)
            coefficient(2) = 2 * (b ** 2.0 - a ** 2.0) * cos(angle) * sin(angle)
            coefficient(3) = (a ** 2.0) * (cos(angle) ** 2.0) + (b ** 2.0) * (sin(angle) ** 2.0)
            coefficient(4) = -2 * coefficient(1) * xc - coefficient(2) * yc
            coefficient(5) = -coefficient(2) * xc - 2 * coefficient(3) * yc
            coefficient(6) = coefficient(1) * (xc ** 2.0) + coefficient(2) * xc * yc + coefficient(3) * (yc ** 2.0) - (a ** 2.0) * (b ** 2.0)

        end subroutine

        !****************************************************************************
        !
        ! Generate a mesh on the surface of an ellipsoid
        ! ellipsoid_info[9]     radius_a    radius_b    radius_c    angle_x   angle_y
        !                       angle_z     center_x_pos   center_y_pos  center_z_pos
        ! longitude_vertex_num (latitide_vertex_num) control the density of the mesh
        ! mesh_info [total_node_num,3]  saves the coordinate of the nodes on the mesh
        !
        !****************************************************************************
        subroutine generate_mesh_on_ellipsoid_surface(ellipsoid_info,longitude_vertex_num,latitude_vertex_num,mesh_info)
            implicit none
            real(kind=8),intent(in)                 ::  ellipsoid_info(9)
            integer(kind=4),intent(in)              ::  longitude_vertex_num, latitude_vertex_num
            real(kind=8),allocatable,intent(out)    ::  mesh_info(:,:)
            integer(kind=4)                         ::  i, j, k, node_num
            real(kind=8)                            ::  longitude_degree, latitude_degree
            real(kind=8)                            ::  trans_rotate_x_axis(4,4), trans_rotate_y_axis(4,4), trans_rotate_z_axis(4,4), trans_rotate(4,4), vertex_pos(4)
            
            if (allocated(mesh_info)) deallocate(mesh_info)
            allocate(mesh_info(longitude_vertex_num*latitude_vertex_num,3))
            mesh_info = 0.0

            !coordinate transform
            trans_rotate_x_axis = 0.0
            trans_rotate_x_axis(1,1) = 1.0
            trans_rotate_x_axis(2,2) = cos(ellipsoid_info(4))
            trans_rotate_x_axis(2,3) = sin(ellipsoid_info(4))
            trans_rotate_x_axis(3,2) = -sin(ellipsoid_info(4))
            trans_rotate_x_axis(3,3) = cos(ellipsoid_info(4))
            trans_rotate_x_axis(4,4) = 1.0
                
            trans_rotate_y_axis = 0.0
            trans_rotate_y_axis(1,1) = cos(ellipsoid_info(5))
            trans_rotate_y_axis(1,3) = -sin(ellipsoid_info(5))
            trans_rotate_y_axis(2,2) = 1
            trans_rotate_y_axis(3,1) = sin(ellipsoid_info(5))
            trans_rotate_y_axis(3,3) = cos(ellipsoid_info(5))
            trans_rotate_y_axis(4,4) = 1.0
                
            trans_rotate_z_axis = 0.0
            trans_rotate_z_axis(1,1) = cos(ellipsoid_info(6))
            trans_rotate_z_axis(1,2) = sin(ellipsoid_info(6))
            trans_rotate_z_axis(2,1) = -sin(ellipsoid_info(6))
            trans_rotate_z_axis(2,2) = cos(ellipsoid_info(6))
            trans_rotate_z_axis(3,3) = 1.0
            trans_rotate_z_axis(4,4) = 1.0
            trans_rotate = matmul(matmul(trans_rotate_x_axis,trans_rotate_y_axis),trans_rotate_z_axis)

            do j = 1, longitude_vertex_num
                longitude_degree = -3.14159 + 2 * 3.14159 / (longitude_vertex_num - 1) * (j - 1)    
                do k = 1, latitude_vertex_num
                    latitude_degree = -3.14159 / 2.0 + 3.14159 / (latitude_vertex_num - 1) * (k - 1)
                    vertex_pos(1) = ellipsoid_info(1) * cos(latitude_degree) * cos(longitude_degree)
                    vertex_pos(2) = ellipsoid_info(2) * cos(latitude_degree) * sin(longitude_degree)
                    vertex_pos(3) = ellipsoid_info(3) * sin(latitude_degree)
                    vertex_pos(4) = 1.0
                    vertex_pos = matmul(vertex_pos,trans_rotate)
                    node_num = (j - 1) * latitude_vertex_num + k
                    mesh_info(node_num,1) = vertex_pos(1)
                    mesh_info(node_num,2) = vertex_pos(2)
                    mesh_info(node_num,3) = vertex_pos(3)
                end do
            end do
        end subroutine
        
        !****************************************************************************
        !
        !  determine whether the mesh on the surface of a geometry is inside a block 
        !  region
        !  mesh_info    [total_node_num,3]  coorinate of all nodes of the mesh
        !
        !****************************************************************************
        function is_geometry_out_of_box_region(mesh_info,left_x_pos,right_x_pos,left_y_pos,right_y_pos,left_z_pos,right_z_pos)
            implicit none
            real(kind=8),allocatable    ::  mesh_info(:,:)
            real(kind=8),intent(in)     ::  left_x_pos, right_x_pos, left_y_pos, right_y_pos, left_z_pos, right_z_pos
            logical                     ::  is_geometry_out_of_box_region
            integer(kind=4)             ::  i, total_node_num

            total_node_num = size(mesh_info,1)

            do i = 1, total_node_num
                if ((mesh_info(i,1) .le. left_x_pos) .or. (mesh_info(i,1) .ge. right_x_pos)) then
                    is_geometry_out_of_box_region = .true.
                    return
                end if
                if ((mesh_info(i,2) .le. left_y_pos) .or. (mesh_info(i,2) .ge. right_y_pos)) then
                    is_geometry_out_of_box_region = .true.
                    return
                end if
                if ((mesh_info(i,3) .le. left_z_pos) .or. (mesh_info(i,3) .ge. right_z_pos)) then
                    is_geometry_out_of_box_region = .true.
                    return
                end if
            end do
            is_geometry_out_of_box_region = .false.
            return
        end function
        
        !****************************************************************************
        !
        !  determine whether the ellipsoid is overlapped with another ellipsoid
        !  detect whether the node of the mesh on the surface of the ellipsoid is
        !  inside another ellipsoid
        !
        !****************************************************************************
        function is_ellipsoid_overlapped(ellipsoid_a_info,ellipsoid_b_info,mesh_a_info,mesh_b_info)
            implicit none
            real(kind=8),intent(in)             ::  ellipsoid_a_info(9), ellipsoid_b_info(9)
            real(kind=8),allocatable,intent(in) ::  mesh_a_info(:,:), mesh_b_info(:,:)
            logical                             ::  is_ellipsoid_overlapped
            integer(kind=4)                     ::  total_node_num_of_mesh_a, total_node_num_of_mesh_b, i
            real(kind=8)                        ::  radius_a, radius_b, radius_c, angle_x, angle_y, angle_z, center_x_pos, center_y_pos, center_z_pos
            real(kind=8)                        ::  ellipsoid_fun_value, trans_rotate_x_axis(4,4), trans_rotate_y_axis(4,4), trans_rotate_z_axis(4,4), trans_rotate(4,4)
            real(kind=8)                        ::  vertex_pos(4)

            total_node_num_of_mesh_a = size(mesh_a_info,1)
            total_node_num_of_mesh_b = size(mesh_b_info,2)

            ! test node of mesh_a inside ellipsoid_b
            radius_a = ellipsoid_b_info(1)
            radius_b = ellipsoid_b_info(2)
            radius_c = ellipsoid_b_info(3)
            angle_x = ellipsoid_b_info(4)
            angle_y = ellipsoid_b_info(5)
            angle_z = ellipsoid_b_info(6)
            center_x_pos = ellipsoid_b_info(7)
            center_y_pos = ellipsoid_b_info(8)
            center_z_pos = ellipsoid_b_info(9)

            trans_rotate_x_axis = 0.0
            trans_rotate_x_axis(1,1) = 1.0
            trans_rotate_x_axis(2,2) = cos(-angle_x)
            trans_rotate_x_axis(2,3) = sin(-angle_x)
            trans_rotate_x_axis(3,2) = -sin(-angle_x)
            trans_rotate_x_axis(3,3) = cos(-angle_x)
            trans_rotate_x_axis(4,4) = 1.0
                
            trans_rotate_y_axis = 0.0
            trans_rotate_y_axis(1,1) = cos(-angle_y)
            trans_rotate_y_axis(1,3) = -sin(-angle_y)
            trans_rotate_y_axis(2,2) = 1
            trans_rotate_y_axis(3,1) = sin(-angle_y)
            trans_rotate_y_axis(3,3) = cos(-angle_y)
            trans_rotate_y_axis(4,4) = 1.0
                
            trans_rotate_z_axis = 0.0
            trans_rotate_z_axis(1,1) = cos(-angle_z)
            trans_rotate_z_axis(1,2) = sin(-angle_z)
            trans_rotate_z_axis(2,1) = -sin(-angle_z)
            trans_rotate_z_axis(2,2) = cos(-angle_z)
            trans_rotate_z_axis(3,3) = 1.0
            trans_rotate_z_axis(4,4) = 1.0
            trans_rotate = matmul(matmul(trans_rotate_z_axis,trans_rotate_y_axis),trans_rotate_x_axis)
            
            do i = 1, total_node_num_of_mesh_a
                vertex_pos(1) = mesh_a_info(i,1) - center_x_pos
                vertex_pos(2) = mesh_a_info(i,2) - center_y_pos
                vertex_pos(3) = mesh_a_info(i,3) - center_z_pos
                vertex_pos(4) = 1.0
                vertex_pos = matmul(vertex_pos,trans_rotate)
                ellipsoid_fun_value = (vertex_pos(1) / radius_a) ** 2.0 + (vertex_pos(2) / radius_b) ** 2.0 + (vertex_pos(3) / radius_c) ** 2.0 - 1.0
                if (ellipsoid_fun_value .lt. 0.0) then
                    is_ellipsoid_overlapped = .true.
                    return
                end if
            end do

            ! test node of mesh_b inside ellipsoid_a
            radius_a = ellipsoid_a_info(1)
            radius_b = ellipsoid_a_info(2)
            radius_c = ellipsoid_a_info(3)
            angle_x = ellipsoid_a_info(4)
            angle_y = ellipsoid_a_info(5)
            angle_z = ellipsoid_a_info(6)
            center_x_pos = ellipsoid_a_info(7)
            center_y_pos = ellipsoid_a_info(8)
            center_z_pos = ellipsoid_a_info(9)

            trans_rotate_x_axis = 0.0
            trans_rotate_x_axis(1,1) = 1.0
            trans_rotate_x_axis(2,2) = cos(-angle_x)
            trans_rotate_x_axis(2,3) = sin(-angle_x)
            trans_rotate_x_axis(3,2) = -sin(-angle_x)
            trans_rotate_x_axis(3,3) = cos(-angle_x)
            trans_rotate_x_axis(4,4) = 1.0
                
            trans_rotate_y_axis = 0.0
            trans_rotate_y_axis(1,1) = cos(-angle_y)
            trans_rotate_y_axis(1,3) = -sin(-angle_y)
            trans_rotate_y_axis(2,2) = 1
            trans_rotate_y_axis(3,1) = sin(-angle_y)
            trans_rotate_y_axis(3,3) = cos(-angle_y)
            trans_rotate_y_axis(4,4) = 1.0
                
            trans_rotate_z_axis = 0.0
            trans_rotate_z_axis(1,1) = cos(-angle_z)
            trans_rotate_z_axis(1,2) = sin(-angle_z)
            trans_rotate_z_axis(2,1) = -sin(-angle_z)
            trans_rotate_z_axis(2,2) = cos(-angle_z)
            trans_rotate_z_axis(3,3) = 1.0
            trans_rotate_z_axis(4,4) = 1.0
            trans_rotate = matmul(matmul(trans_rotate_z_axis,trans_rotate_y_axis),trans_rotate_x_axis)
            
            do i = 1, total_node_num_of_mesh_b
                vertex_pos(1) = mesh_b_info(i,1) - center_x_pos
                vertex_pos(2) = mesh_b_info(i,2) - center_y_pos
                vertex_pos(3) = mesh_b_info(i,3) - center_z_pos
                vertex_pos(4) = 1.0
                vertex_pos = matmul(vertex_pos,trans_rotate)
                ellipsoid_fun_value = (vertex_pos(1) / radius_a) ** 2.0 + (vertex_pos(2) / radius_b) ** 2.0 + (vertex_pos(3) / radius_c) ** 2.0 - 1.0
                if (ellipsoid_fun_value .lt. 0.0) then
                    is_ellipsoid_overlapped = .true.
                    return
                end if
            end do
            
            is_ellipsoid_overlapped = .false.
            return
        end function

        !****************************************************************************
        !
        ! Generate a mesh on the surface of a cylinder
        ! cylinder_info[9]     radius_a    radius_b    height    angle_x   angle_y
        !                       angle_z     center_x_pos   center_y_pos  center_z_pos
        ! total_angle_segment (total_height_segment) control the density of the mesh
        ! mesh_info [total_node_num,3]  saves the coordinate of the nodes on the mesh
        !
        !****************************************************************************
        subroutine generate_mesh_on_cylinder_surface(cylinder_info,total_angle_segment,total_height_segment,mesh_info)
            implicit none
            real(kind=8),intent(in)                 ::  cylinder_info(9)
            integer(kind=4),intent(in)              ::  total_angle_segment, total_height_segment
            real(kind=8),allocatable,intent(out)    ::  mesh_info(:,:)
            integer(kind=4)                         ::  i, j, k, node_num
            real(kind=8)                            ::  angle, point_pos(4)
            real(kind=8)                            ::  trans_rotate_x_axis(4,4), trans_rotate_y_axis(4,4), trans_rotate_z_axis(4,4), trans_rotate(4,4)
            
            if (allocated(mesh_info)) deallocate(mesh_info)
            allocate(mesh_info(total_angle_segment*total_height_segment,3))
            mesh_info = 0.0

            ! rotation matrix
            trans_rotate_x_axis = 0.0
            trans_rotate_x_axis(1,1) =  1.0
            trans_rotate_x_axis(2,2) =  cos(cylinder_info(4))
            trans_rotate_x_axis(2,3) =  sin(cylinder_info(4))
            trans_rotate_x_axis(3,2) = -sin(cylinder_info(4))
            trans_rotate_x_axis(3,3) =  cos(cylinder_info(4))
            trans_rotate_x_axis(4,4) =  1.0
                    
            trans_rotate_y_axis = 0.0
            trans_rotate_y_axis(1,1) =  cos(cylinder_info(5))
            trans_rotate_y_axis(1,3) = -sin(cylinder_info(5))
            trans_rotate_y_axis(2,2) =  1.0
            trans_rotate_y_axis(3,1) =  sin(cylinder_info(5))
            trans_rotate_y_axis(3,3) =  cos(cylinder_info(5))
            trans_rotate_y_axis(4,4) =  1.0
                    
            trans_rotate_z_axis = 0.0
            trans_rotate_z_axis(1,1) =  cos(cylinder_info(6))
            trans_rotate_z_axis(1,2) =  sin(cylinder_info(6))
            trans_rotate_z_axis(2,1) = -sin(cylinder_info(6))
            trans_rotate_z_axis(2,2) =  cos(cylinder_info(6))
            trans_rotate_z_axis(3,3) =  1.0
            trans_rotate_z_axis(4,4) =  1.0
            trans_rotate = matmul(matmul(trans_rotate_x_axis,trans_rotate_y_axis),trans_rotate_z_axis)

            do i = 1, total_angle_segment
                angle = -3.14159 + 2.0 * 3.14159 / (total_angle_segment - 1) * (i - 1)
                do j = 1, total_height_segment   
                    point_pos(1) = cylinder_info(1) * cos(angle)
                    point_pos(2) = cylinder_info(2) * sin(angle)
                    point_pos(3) = -0.5 * cylinder_info(3) + (j - 1) * cylinder_info(3) / (total_height_segment - 1)
                    point_pos(4) = 1.0
                    point_pos = matmul(point_pos,trans_rotate)

                    node_num = (i - 1) * total_height_segment + j
                    mesh_info(node_num,1) = point_pos(1) + cylinder_info(7)
                    mesh_info(node_num,2) = point_pos(2) + cylinder_info(8)
                    mesh_info(node_num,3) = point_pos(3) + cylinder_info(9)
           
                end do
            end do

        end subroutine

        !****************************************************************************
        !
        ! Generate a mesh on the surface of a capsule
        ! capsule_info[8]       radius    height    angle_x   angle_y     angle_z
        !                       center_x_pos   center_y_pos  center_z_pos
        ! total_height_segment (total_alpha_segment,total_beta_segment) 
        ! control the density of the mesh
        ! mesh_info [total_node_num,3]  saves the coordinate of the nodes on the mesh
        !
        !****************************************************************************
        subroutine generate_mesh_on_capsule_surface(capsule_info,total_alpha_segment,total_beta_segment,total_height_segment,mesh_info)
            implicit none
            real(kind=8),intent(in)                 ::  capsule_info(8)
            integer(kind=4),intent(in)              ::  total_alpha_segment, total_beta_segment, total_height_segment
            real(kind=8),allocatable,intent(out)    ::  mesh_info(:,:)
            integer(kind=4)                         ::  i, j, node_num, total_node_num
            real(kind=8)                            ::  alpha, beta, radius, point_pos(4)
            real(kind=8)                            ::  trans_rotate_x_axis(4,4), trans_rotate_y_axis(4,4), trans_rotate_z_axis(4,4), trans_rotate(4,4)
            
            total_node_num = (total_height_segment - 2) * total_alpha_segment + 2 * total_alpha_segment * total_beta_segment
            if (allocated(mesh_info)) deallocate(mesh_info)
            allocate(mesh_info(total_node_num,3))
            mesh_info = 0.0

            ! rotation matrix
            trans_rotate_x_axis = 0.0
            trans_rotate_x_axis(1,1) =  1.0
            trans_rotate_x_axis(2,2) =  cos(capsule_info(3))
            trans_rotate_x_axis(2,3) =  sin(capsule_info(3))
            trans_rotate_x_axis(3,2) = -sin(capsule_info(3))
            trans_rotate_x_axis(3,3) =  cos(capsule_info(3))
            trans_rotate_x_axis(4,4) =  1.0
                    
            trans_rotate_y_axis = 0.0
            trans_rotate_y_axis(1,1) =  cos(capsule_info(4))
            trans_rotate_y_axis(1,3) = -sin(capsule_info(4))
            trans_rotate_y_axis(2,2) =  1.0
            trans_rotate_y_axis(3,1) =  sin(capsule_info(4))
            trans_rotate_y_axis(3,3) =  cos(capsule_info(4))
            trans_rotate_y_axis(4,4) =  1.0
                    
            trans_rotate_z_axis = 0.0
            trans_rotate_z_axis(1,1) =  cos(capsule_info(5))
            trans_rotate_z_axis(1,2) =  sin(capsule_info(5))
            trans_rotate_z_axis(2,1) = -sin(capsule_info(5))
            trans_rotate_z_axis(2,2) =  cos(capsule_info(5))
            trans_rotate_z_axis(3,3) =  1.0
            trans_rotate_z_axis(4,4) =  1.0
            trans_rotate = matmul(matmul(trans_rotate_x_axis,trans_rotate_y_axis),trans_rotate_z_axis)

            ! spherical tips of capsule
            ! bottom
            do i = 1, total_alpha_segment
                alpha = -3.14159 + 2 * 3.14159 / (total_alpha_segment - 1) * (i - 1)    
                do j = 1, total_beta_segment
                    beta = -3.14159 / 2.0 + 3.14159 / 2.0 / (total_beta_segment - 1) * (j - 1)
                    point_pos(1) = capsule_info(1) * cos(beta) * cos(alpha)
                    point_pos(2) = capsule_info(1) * cos(beta) * sin(alpha)
                    point_pos(3) = capsule_info(1) * sin(beta) - 0.5 * capsule_info(2)
                    point_pos(4) = 1.0
                    point_pos = matmul(point_pos,trans_rotate)
                    node_num = (i - 1) * total_beta_segment + j
                    mesh_info(node_num,1) = point_pos(1) + capsule_info(6)
                    mesh_info(node_num,2) = point_pos(2) + capsule_info(7)
                    mesh_info(node_num,3) = point_pos(3) + capsule_info(8)
                end do
            end do
            ! top
            do i = 1, total_alpha_segment
                alpha = -3.14159 + 2 * 3.14159 / (total_alpha_segment - 1) * (i - 1)    
                do j = 1, total_beta_segment
                    beta = 3.14159 / 2.0 / (total_beta_segment - 1) * (j - 1)
                    point_pos(1) = capsule_info(1) * cos(beta) * cos(alpha)
                    point_pos(2) = capsule_info(1) * cos(beta) * sin(alpha)
                    point_pos(3) = capsule_info(1) * sin(beta) + 0.5 * capsule_info(2)
                    point_pos(4) = 1.0
                    point_pos = matmul(point_pos,trans_rotate)
                    node_num = total_alpha_segment * total_beta_segment + (i - 1) * total_beta_segment + j
                    mesh_info(node_num,1) = point_pos(1) + capsule_info(6)
                    mesh_info(node_num,2) = point_pos(2) + capsule_info(7)
                    mesh_info(node_num,3) = point_pos(3) + capsule_info(8)
                end do
            end do

            ! cylinder part of capsule
            do i = 1, total_alpha_segment
                alpha = -3.14159 + 2 * 3.14159 / (total_alpha_segment - 1) * (i - 1)
                do j = 2, total_height_segment - 1              
                    point_pos(1) = capsule_info(1) * cos(alpha)
                    point_pos(2) = capsule_info(1) * sin(alpha)
                    point_pos(3) = -0.5 * capsule_info(2) + (j - 1) * capsule_info(2) / (total_height_segment - 1)
                    point_pos(4) = 1.0
                    point_pos = matmul(point_pos,trans_rotate)
                    node_num = 2 * total_alpha_segment * total_beta_segment + (i - 1) * (total_height_segment - 2) + j - 1
                    mesh_info(node_num,1) = point_pos(1) + capsule_info(6)
                    mesh_info(node_num,2) = point_pos(2) + capsule_info(7)
                    mesh_info(node_num,3) = point_pos(3) + capsule_info(8)   
                end do
            end do

        end subroutine

        !****************************************************************************
        !
        !  determine whether the cylinder is overlapped with another cylinder
        !  detect whether the node of the mesh on the surface of the cylinder is
        !  inside another cylinder
        !
        !****************************************************************************
        function is_cylinder_overlapped(cylinder_a_info,cylinder_b_info,mesh_a_info,mesh_b_info)
            implicit none
            real(kind=8),intent(in)             ::  cylinder_a_info(9), cylinder_b_info(9)
            real(kind=8),allocatable,intent(in) ::  mesh_a_info(:,:), mesh_b_info(:,:)
            logical                             ::  is_cylinder_overlapped
            integer(kind=4)                     ::  total_node_num_of_mesh_a, total_node_num_of_mesh_b, i
            real(kind=8)                        ::  radius_a, radius_b, height, angle_x, angle_y, angle_z, center_x_pos, center_y_pos, center_z_pos, ellipse_fun_value
            real(kind=8)                        ::  cylinder_fun_value, trans_rotate_x_axis(4,4), trans_rotate_y_axis(4,4), trans_rotate_z_axis(4,4), trans_rotate(4,4)
            real(kind=8)                        ::  point_pos(4)

            total_node_num_of_mesh_a = size(mesh_a_info,1)
            total_node_num_of_mesh_b = size(mesh_b_info,2)

            ! test node of mesh_a inside cylinder_b
            radius_a = cylinder_b_info(1)
            radius_b = cylinder_b_info(2)
            height = cylinder_b_info(3)
            angle_x = cylinder_b_info(4)
            angle_y = cylinder_b_info(5)
            angle_z = cylinder_b_info(6)
            center_x_pos = cylinder_b_info(7)
            center_y_pos = cylinder_b_info(8)
            center_z_pos = cylinder_b_info(9)

            trans_rotate_x_axis = 0.0
            trans_rotate_x_axis(1,1) = 1.0
            trans_rotate_x_axis(2,2) = cos(-angle_x)
            trans_rotate_x_axis(2,3) = sin(-angle_x)
            trans_rotate_x_axis(3,2) = -sin(-angle_x)
            trans_rotate_x_axis(3,3) = cos(-angle_x)
            trans_rotate_x_axis(4,4) = 1.0
                
            trans_rotate_y_axis = 0.0
            trans_rotate_y_axis(1,1) = cos(-angle_y)
            trans_rotate_y_axis(1,3) = -sin(-angle_y)
            trans_rotate_y_axis(2,2) = 1
            trans_rotate_y_axis(3,1) = sin(-angle_y)
            trans_rotate_y_axis(3,3) = cos(-angle_y)
            trans_rotate_y_axis(4,4) = 1.0
                
            trans_rotate_z_axis = 0.0
            trans_rotate_z_axis(1,1) = cos(-angle_z)
            trans_rotate_z_axis(1,2) = sin(-angle_z)
            trans_rotate_z_axis(2,1) = -sin(-angle_z)
            trans_rotate_z_axis(2,2) = cos(-angle_z)
            trans_rotate_z_axis(3,3) = 1.0
            trans_rotate_z_axis(4,4) = 1.0
            trans_rotate = matmul(matmul(trans_rotate_z_axis,trans_rotate_y_axis),trans_rotate_x_axis)
            
            do i = 1, total_node_num_of_mesh_a
                point_pos(1) = mesh_a_info(i,1) - center_x_pos
                point_pos(2) = mesh_a_info(i,2) - center_y_pos
                point_pos(3) = mesh_a_info(i,3) - center_z_pos
                point_pos(4) = 1.0
                point_pos = matmul(point_pos,trans_rotate)
                if ((point_pos(3) .ge. -cylinder_b_info(3) * 0.5) .and. (point_pos(3) .le. cylinder_b_info(3) * 0.5)) then
                    ellipse_fun_value = (point_pos(1) / radius_a) ** 2.0 + (point_pos(2) / radius_b) ** 2.0 - 1.0
                    if (ellipse_fun_value .le. 0.0) then
                        is_cylinder_overlapped = .true.
                        return
                    end if
                end if    
            end do

            ! test node of mesh_b inside cylinder_a
            radius_a = cylinder_a_info(1)
            radius_b = cylinder_a_info(2)
            height = cylinder_a_info(3)
            angle_x = cylinder_a_info(4)
            angle_y = cylinder_a_info(5)
            angle_z = cylinder_a_info(6)
            center_x_pos = cylinder_a_info(7)
            center_y_pos = cylinder_a_info(8)
            center_z_pos = cylinder_a_info(9)

            trans_rotate_x_axis = 0.0
            trans_rotate_x_axis(1,1) = 1.0
            trans_rotate_x_axis(2,2) = cos(-angle_x)
            trans_rotate_x_axis(2,3) = sin(-angle_x)
            trans_rotate_x_axis(3,2) = -sin(-angle_x)
            trans_rotate_x_axis(3,3) = cos(-angle_x)
            trans_rotate_x_axis(4,4) = 1.0
                
            trans_rotate_y_axis = 0.0
            trans_rotate_y_axis(1,1) = cos(-angle_y)
            trans_rotate_y_axis(1,3) = -sin(-angle_y)
            trans_rotate_y_axis(2,2) = 1
            trans_rotate_y_axis(3,1) = sin(-angle_y)
            trans_rotate_y_axis(3,3) = cos(-angle_y)
            trans_rotate_y_axis(4,4) = 1.0
                
            trans_rotate_z_axis = 0.0
            trans_rotate_z_axis(1,1) = cos(-angle_z)
            trans_rotate_z_axis(1,2) = sin(-angle_z)
            trans_rotate_z_axis(2,1) = -sin(-angle_z)
            trans_rotate_z_axis(2,2) = cos(-angle_z)
            trans_rotate_z_axis(3,3) = 1.0
            trans_rotate_z_axis(4,4) = 1.0
            trans_rotate = matmul(matmul(trans_rotate_z_axis,trans_rotate_y_axis),trans_rotate_x_axis)
            
            do i = 1, total_node_num_of_mesh_b
                point_pos(1) = mesh_b_info(i,1) - center_x_pos
                point_pos(2) = mesh_b_info(i,2) - center_y_pos
                point_pos(3) = mesh_b_info(i,3) - center_z_pos
                point_pos(4) = 1.0
                point_pos = matmul(point_pos,trans_rotate)
                if ((point_pos(3) .ge. -cylinder_a_info(3) * 0.5) .and. (point_pos(3) .le. cylinder_a_info(3) * 0.5)) then
                    ellipse_fun_value = (point_pos(1) / radius_a) ** 2.0 + (point_pos(2) / radius_b) ** 2.0 - 1.0
                    if (ellipse_fun_value .le. 0.0) then
                        is_cylinder_overlapped = .true.
                        return
                    end if
                end if    
            end do
            
            is_cylinder_overlapped = .false.
            return
        end function

        !****************************************************************************
        !
        !  determine whether the capsule is overlapped with another capsule
        !  detect whether the node of the mesh on the surface of the capsule is
        !  inside another capsule
        !
        !****************************************************************************
        function is_capsule_overlapped(capsule_a_info,capsule_b_info,mesh_a_info,mesh_b_info)
            implicit none
            real(kind=8),intent(in)             ::  capsule_a_info(8), capsule_b_info(8)
            real(kind=8),allocatable,intent(in) ::  mesh_a_info(:,:), mesh_b_info(:,:)
            logical                             ::  is_capsule_overlapped
            integer(kind=4)                     ::  total_node_num_of_mesh_a, total_node_num_of_mesh_b, i
            real(kind=8)                        ::  radius, height, angle_x, angle_y, angle_z, center_x_pos, center_y_pos, center_z_pos, value
            real(kind=8)                        ::  capsule_fun_value, trans_rotate_x_axis(4,4), trans_rotate_y_axis(4,4), trans_rotate_z_axis(4,4), trans_rotate(4,4)
            real(kind=8)                        ::  point_pos(4)

            total_node_num_of_mesh_a = size(mesh_a_info,1)
            total_node_num_of_mesh_b = size(mesh_b_info,2)

            ! test node of mesh_a inside capsule_b
            radius = capsule_b_info(1)
            height = capsule_b_info(2)
            angle_x = capsule_b_info(3)
            angle_y = capsule_b_info(4)
            angle_z = capsule_b_info(5)
            center_x_pos = capsule_b_info(6)
            center_y_pos = capsule_b_info(7)
            center_z_pos = capsule_b_info(8)

            trans_rotate_x_axis = 0.0
            trans_rotate_x_axis(1,1) = 1.0
            trans_rotate_x_axis(2,2) = cos(-angle_x)
            trans_rotate_x_axis(2,3) = sin(-angle_x)
            trans_rotate_x_axis(3,2) = -sin(-angle_x)
            trans_rotate_x_axis(3,3) = cos(-angle_x)
            trans_rotate_x_axis(4,4) = 1.0
                
            trans_rotate_y_axis = 0.0
            trans_rotate_y_axis(1,1) = cos(-angle_y)
            trans_rotate_y_axis(1,3) = -sin(-angle_y)
            trans_rotate_y_axis(2,2) = 1
            trans_rotate_y_axis(3,1) = sin(-angle_y)
            trans_rotate_y_axis(3,3) = cos(-angle_y)
            trans_rotate_y_axis(4,4) = 1.0
                
            trans_rotate_z_axis = 0.0
            trans_rotate_z_axis(1,1) = cos(-angle_z)
            trans_rotate_z_axis(1,2) = sin(-angle_z)
            trans_rotate_z_axis(2,1) = -sin(-angle_z)
            trans_rotate_z_axis(2,2) = cos(-angle_z)
            trans_rotate_z_axis(3,3) = 1.0
            trans_rotate_z_axis(4,4) = 1.0
            trans_rotate = matmul(matmul(trans_rotate_z_axis,trans_rotate_y_axis),trans_rotate_x_axis)
            
            do i = 1, total_node_num_of_mesh_a
                point_pos(1) = mesh_a_info(i,1) - center_x_pos
                point_pos(2) = mesh_a_info(i,2) - center_y_pos
                point_pos(3) = mesh_a_info(i,3) - center_z_pos
                point_pos(4) = 1.0
                point_pos = matmul(point_pos,trans_rotate)

                ! whether in the cylinder part of the capsule
                if ((point_pos(3) .ge. -capsule_b_info(2) * 0.5) .and. (point_pos(3) .le. capsule_b_info(2) * 0.5)) then
                    value = (point_pos(1) ** 2.0 + point_pos(2) ** 2.0) - radius ** 2.0
                    if (value .le. 0.0) then
                        is_capsule_overlapped = .true.
                        return
                    end if
                end if

                ! whether in the spherical tips of the capsule
                value = (point_pos(1) ** 2.0 + point_pos(2) ** 2.0 + (point_pos(3) - (-0.5 * capsule_b_info(2))) ** 2.0) - radius ** 2.0
                if (value .lt. 0.0) then
                    is_capsule_overlapped = .true.
                    return
                end if
                value = (point_pos(1) ** 2.0 + point_pos(2) ** 2.0 + (point_pos(3) - (0.5 * capsule_b_info(2))) ** 2.0) - radius ** 2.0
                if (value .lt. 0.0) then
                    is_capsule_overlapped = .true.
                    return
                end if
            end do

            ! test node of mesh_b inside capsule_a
            radius = capsule_a_info(1)
            height = capsule_a_info(2)
            angle_x = capsule_a_info(3)
            angle_y = capsule_a_info(4)
            angle_z = capsule_a_info(5)
            center_x_pos = capsule_a_info(6)
            center_y_pos = capsule_a_info(7)
            center_z_pos = capsule_a_info(8)

            trans_rotate_x_axis = 0.0
            trans_rotate_x_axis(1,1) = 1.0
            trans_rotate_x_axis(2,2) = cos(-angle_x)
            trans_rotate_x_axis(2,3) = sin(-angle_x)
            trans_rotate_x_axis(3,2) = -sin(-angle_x)
            trans_rotate_x_axis(3,3) = cos(-angle_x)
            trans_rotate_x_axis(4,4) = 1.0
                
            trans_rotate_y_axis = 0.0
            trans_rotate_y_axis(1,1) = cos(-angle_y)
            trans_rotate_y_axis(1,3) = -sin(-angle_y)
            trans_rotate_y_axis(2,2) = 1
            trans_rotate_y_axis(3,1) = sin(-angle_y)
            trans_rotate_y_axis(3,3) = cos(-angle_y)
            trans_rotate_y_axis(4,4) = 1.0
                
            trans_rotate_z_axis = 0.0
            trans_rotate_z_axis(1,1) = cos(-angle_z)
            trans_rotate_z_axis(1,2) = sin(-angle_z)
            trans_rotate_z_axis(2,1) = -sin(-angle_z)
            trans_rotate_z_axis(2,2) = cos(-angle_z)
            trans_rotate_z_axis(3,3) = 1.0
            trans_rotate_z_axis(4,4) = 1.0
            trans_rotate = matmul(matmul(trans_rotate_z_axis,trans_rotate_y_axis),trans_rotate_x_axis)
            
            do i = 1, total_node_num_of_mesh_b
                point_pos(1) = mesh_b_info(i,1) - center_x_pos
                point_pos(2) = mesh_b_info(i,2) - center_y_pos
                point_pos(3) = mesh_b_info(i,3) - center_z_pos
                point_pos(4) = 1.0
                point_pos = matmul(point_pos,trans_rotate)

                ! whether in the cylinder part of the capsule
                if ((point_pos(3) .ge. -capsule_a_info(2) * 0.5) .and. (point_pos(3) .le. capsule_a_info(2) * 0.5)) then
                    value = (point_pos(1) ** 2.0 + point_pos(2) ** 2.0) - radius ** 2.0
                    if (value .le. 0.0) then
                        is_capsule_overlapped = .true.
                        return
                    end if
                end if

                ! whether in the spherical tips of the capsule
                value = (point_pos(1) ** 2.0 + point_pos(2) ** 2.0 + (point_pos(3) - (-0.5 * capsule_a_info(2))) ** 2.0) - radius ** 2.0
                if (value .lt. 0.0) then
                    is_capsule_overlapped = .true.
                    return
                end if
                value = (point_pos(1) ** 2.0 + point_pos(2) ** 2.0 + (point_pos(3) - (0.5 * capsule_a_info(2))) ** 2.0) - radius ** 2.0
                if (value .lt. 0.0) then
                    is_capsule_overlapped = .true.
                    return
                end if  
            end do
            
            is_capsule_overlapped = .false.
            return
        end function

end module

!dec$ endif