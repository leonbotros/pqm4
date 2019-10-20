.syntax unified
.cpu cortex-m4
.thumb

.global asm_gf2x_mul1
.type asm_gf2x_mul1, %function
.align 2
asm_gf2x_mul1:
push.w {r4-r7, lr}
mov r7, #1
and r6, r1, r7
mul r4, r6, r2
mov r3, #0
and r6, r1, r7, lsl #1
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #2
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #3
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #4
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #5
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #6
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #7
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #8
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #9
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #10
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #11
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #12
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #13
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #14
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #15
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #16
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #17
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #18
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #19
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #20
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #21
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #22
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #23
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #24
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #25
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #26
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #27
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #28
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #29
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #30
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
and r6, r1, r7, lsl #31
umull r6, r5, r6, r2
eor r4, r6
eor r3, r5
stm r0, {r3, r4}
pop {r4-r7, pc}
