.syntax unified
.cpu cortex-m4
.thumb

.global asm_gf2x_shift_right
.type asm_gf2x_shift_right, %function
.align 2
asm_gf2x_shift_right:
    push {r4-r11, lr}

    /* set zero carry */
    adds.w r1, #0

    /* shift first 5 words */
    ldm r0, {r1-r5}
    rrxs r1, r1
    rrxs r2, r2
    rrxs r3, r3
    rrxs r4, r4
    rrxs r5, r5
    stm.w r0!, {r1-r5}

    /* shift the remainder of words */
    b.w shift_right


.global asm_gf2x_rotate_right
.type asm_gf2x_rotate_right, %function
.align 2
asm_gf2x_rotate_right:
    push {r4-r11, lr}

    /* set lsb as msb and shift 5 words*/
    movw r8, #6516
    ldr.w r9, [r0, r8]
    ldm r0, {r1-r5}
    and r9, r9, #1
    adds.w r1, #0
    rrxs r1, r1
    orr r1, r1, r9, lsl #18 /* msb is at index 18 in most significant word */
    rrxs r2, r2
    rrxs r3, r3
    rrxs r4, r4
    rrxs r5, r5
    stm.w r0!, {r1-r5}

    /* shift remainder of words */
    b.w shift_right

shift_right:
    .rept 125
        ldm r0, {r1-r12,r14}
        rrxs r1, r1
        rrxs r2, r2
        rrxs r3, r3
        rrxs r4, r4
        rrxs r5, r5
        rrxs r6, r6
        rrxs r7, r7
        rrxs r8, r8
        rrxs r9, r9
        rrxs r10, r10
        rrxs r11, r11
        rrxs r12, r12
        rrxs r14, r14
        stm r0!, {r1-r12,r14}
    .endr

    b.w shift_end

shift_end:
    pop {r4-r11, pc}


/* asm_gf2x_cswap(DIGIT *a, DIGIT *b, int mask) */
.global asm_gf2x_cswap
.type asm_gf2x_cswap, %function
.align 2
asm_gf2x_cswap:
    push {r4-r11, lr}
    mov r14, #543
    mov.w r12, r2

    loop1:
        ldm r0, {r2-r4}
        ldm r1, {r5-r7}

        /* t = mask & (a ^ b) */
        eor r8, r2, r5
        eor r9, r3, r6
        eor r10, r4, r7

        and r8, r8, r12
        and r9, r9, r12
        and r10, r10, r12

        /* a = a ^ t */
        eor r2, r2, r8
        eor r3, r3, r9
        eor r4, r4, r10

        /* b = b ^ t */
        eor r5, r5, r8
        eor r6, r6, r9
        eor r7, r7, r10

        stm.w r0!, {r2-r4}
        stm.w r1!, {r5-r7}

        subs.w r14, #1
        bne.w loop1

    ldr.w r2, [r0]
    ldr.w r3, [r1]
    eor r4, r2, r3
    and r4, r4, r12
    eor r2, r2, r4
    eor r3, r3, r4
    str.w r2, [r0]
    str.w r3, [r1]

    pop {r4-r11, pc}


.global asm_gf2x_mult_acc_scalar
.type asm_gf2x_mult_acc_scalar, %function
.align 2
asm_gf2x_mult_acc_scalar:
    push {r4-r11, lr}
    mov r14, #326

    loop2:
        ldm r0, {r3-r7}
        ldm r1!, {r8-r12}

        /* poly1 = poly1 ^ (mask & poly2) */
        and r8, r8, r2
        and r9, r9, r2
        and r10, r10, r2
        and r11, r11, r2
        and r12, r12, r2

        eor r3, r3, r8
        eor r4, r4, r9
        eor r5, r5, r10
        eor r6, r6, r11
        eor r7, r7, r12

        stm.w r0!, {r3-r7}
        subs.w r14, #1
        bne.w loop2

    pop {r4-r11, pc}


/* void gf2x_add(DIGIT *R, const DIGIT *A, const DIGIT *B, size_t len) */
.global asm_gf2x_add
.type asm_gf2x_add, %function
.align 2
asm_gf2x_add:
    push {r4-r11, lr}
    mov.w r14, r3

    loop3:
        cmp r14, #5
        blt remlt5
        ldm r1!, {r3-r7}
        ldm r2!, {r8-r12}
        eor r3, r8
        eor r4, r9
        eor r5, r10
        eor r6, r11
        eor r7, r12
        stm.w r0!, {r3-r7}
        subs r14, #5
        b.w loop3
    
    remlt5:
        cmp r14, #4
        beq addrem4
        cmp r14, #3
        beq addrem3
        cmp r14, #2
        beq addrem2
        cmp r14, #1
        beq addrem1
        b add_end

    addrem4:
        ldm r1, {r3-r6}
        ldm r2, {r7-r10}
        eor r3, r7
        eor r4, r8
        eor r5, r9
        eor r6, r10
        stm r0, {r3-r6}
        b add_end

    addrem3:
        ldm r1, {r3-r5}
        ldm r2, {r6-r8}
        eor r3, r6
        eor r4, r7
        eor r5, r8
        stm r0, {r3-r5}
        b add_end

    addrem2:
        ldm r1, {r3-r4}
        ldm r2, {r5-r6}
        eor r3, r5
        eor r4, r6
        stm r0, {r3-r4}
        b add_end

    addrem1:
        ldr r3, [r1]
        ldr r4, [r2]
        eor r3, r4
        str r3, [r0]
        b add_end

    add_end:
        pop {r4-r11, pc}
