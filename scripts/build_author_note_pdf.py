#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.platypus import (
    BaseDocTemplate,
    Frame,
    PageBreak,
    PageTemplate,
    Paragraph,
    Spacer,
)


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "output" / "pdf" / "treeeval_mv_parameter_note.pdf"


def footer(canvas, document) -> None:
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#C9CED6"))
    canvas.setLineWidth(0.4)
    canvas.line(0.62 * inch, 0.47 * inch, 7.88 * inch, 0.47 * inch)
    canvas.setFillColor(colors.HexColor("#596273"))
    canvas.setFont("Helvetica", 7.5)
    canvas.drawString(0.62 * inch, 0.29 * inch, "Draft for author feedback")
    canvas.drawRightString(7.88 * inch, 0.29 * inch, f"{document.page}")
    canvas.restoreState()


def build() -> None:
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    doc = BaseDocTemplate(
        str(OUTPUT),
        pagesize=letter,
        leftMargin=0.62 * inch,
        rightMargin=0.62 * inch,
        topMargin=0.55 * inch,
        bottomMargin=0.58 * inch,
        title="Questions on the matching-vector TreeEval construction",
        author="Max",
    )
    frame = Frame(
        doc.leftMargin,
        doc.bottomMargin,
        doc.width,
        doc.height,
        leftPadding=0,
        rightPadding=0,
        topPadding=0,
        bottomPadding=0,
    )
    doc.addPageTemplates(PageTemplate(id="note", frames=[frame], onPage=footer))

    styles = getSampleStyleSheet()
    title = ParagraphStyle(
        "NoteTitle",
        parent=styles["Title"],
        fontName="Helvetica-Bold",
        fontSize=16.5,
        leading=19,
        alignment=TA_CENTER,
        textColor=colors.HexColor("#172033"),
        spaceAfter=5,
    )
    subtitle = ParagraphStyle(
        "Subtitle",
        parent=styles["Normal"],
        fontName="Helvetica",
        fontSize=8.7,
        leading=11,
        alignment=TA_CENTER,
        textColor=colors.HexColor("#596273"),
        spaceAfter=10,
    )
    heading = ParagraphStyle(
        "Heading",
        parent=styles["Heading2"],
        fontName="Helvetica-Bold",
        fontSize=11.2,
        leading=13.5,
        textColor=colors.HexColor("#1F4E79"),
        spaceBefore=7,
        spaceAfter=4,
        keepWithNext=True,
    )
    body = ParagraphStyle(
        "Body",
        parent=styles["BodyText"],
        fontName="Times-Roman",
        fontSize=9.15,
        leading=11.35,
        alignment=TA_LEFT,
        textColor=colors.HexColor("#20242B"),
        spaceAfter=4.2,
    )
    formula = ParagraphStyle(
        "Formula",
        parent=body,
        fontName="Courier",
        fontSize=8.4,
        leading=10.4,
        leftIndent=15,
        rightIndent=15,
        alignment=TA_CENTER,
        backColor=colors.HexColor("#F3F6FA"),
        borderColor=colors.HexColor("#D7DEE8"),
        borderWidth=0.4,
        borderPadding=4,
        spaceBefore=2,
        spaceAfter=5,
    )
    refs = ParagraphStyle(
        "Refs",
        parent=body,
        fontSize=7.7,
        leading=9.2,
        leftIndent=12,
        firstLineIndent=-12,
        spaceAfter=2,
    )

    story = [
        Paragraph("Questions on the matching-vector TreeEval construction", title),
        Paragraph("Notes for the authors | 1 August 2026", subtitle),
        Paragraph(
            "I wrote exact checkers for the matching-vector and selector identities while working through the paper. The main question that came out of this is about the binary representation of the catalytic registers. I have included two smaller parameter observations for context.",
            body,
        ),
        Paragraph("1. Can the catalytic registers be packed jointly?", heading),
        Paragraph(
            "HPR Remark 2.3 validates <font name='Courier'>Z_m</font> coordinates separately and gives <font name='Courier'>O(d log(dm))</font> binary catalytic space. It seems possible to avoid the extra <font name='Courier'>log d</font> factor by packing all <font name='Courier'>D=Theta(d)</font> coordinates into one integer.",
            body,
        ),
        Paragraph("M=m^D,  b=ceil(log_2 M),  y=qM+x,  q in {0,1},  0&lt;=x&lt;M", formula),
        Paragraph(
            "Interpret the arbitrary initial <font name='Courier'>b</font> catalyst bits as <font name='Courier'>y in [0,2^b)</font>. Since <font name='Courier'>2^b&lt;2M</font>, save the one-bit quotient <font name='Courier'>q</font>, subtract <font name='Courier'>qM</font> in place, and read <font name='Courier'>x=sum_i x_i m^i</font> as the vector in <font name='Courier'>Z_m^D</font>. Once the ring algorithm has restored <font name='Courier'>x</font>, adding <font name='Courier'>qM</font> restores the original binary tape exactly.",
            body,
        ),
        Paragraph("Packed-coordinate access", heading),
        Paragraph(
            "Coordinate <font name='Courier'>i</font> is <font name='Courier'>floor(x/m^i) mod m</font>. If it is changed by <font name='Courier'>a</font>, retain the old <font name='Courier'>O(log m)</font>-bit digit and add the signed integer",
            body,
        ),
        Paragraph("((x_i+a mod m)-x_i)m^i", formula),
        Paragraph(
            "to the packed tape in a low-to-high carry or borrow pass. Bits of powers, products, and quotients can be generated in logspace using the uniform TC0 arithmetic of Hesse-Allender-Barrington. Reads and inner products stream over the packed digits. HPR use a constant number of registers, so their coordinates can be concatenated and register swaps handled by a constant-size logical permutation.",
            body,
        ),
        Paragraph(
            "This appears to use <font name='Courier'>ceil(D log_2 m)</font> catalyst bits and <font name='Courier'>O(log m+log(D log m))</font> additional work bits, with polynomial overhead. In Theorem 3.1 it would replace <font name='Courier'>O(d log(dm))</font> by <font name='Courier'>O(d log m)</font> without changing the stated free-space or polynomial-time bounds.",
            body,
        ),
        Paragraph(
            "BCKLS remark after Lemma 15 that stronger compression of the high-order bits should reduce their non-power-of-two simulation to linear catalyst size. I do not view the target as new; the question is whether this state map is a valid realization for HPR's product-ring registers, or whether packed coordinate access violates a constraint of the machine model.",
            body,
        ),
        Paragraph("Checks", heading),
        Paragraph(
            "I checked the state map in an integer implementation and a separate destructive bit-tape implementation. For <font name='Courier'>m=3,D=4</font>, both exhaust all 128 initial seven-bit tapes, every coordinate, and every modular update: 1,536 update/inverse transitions in each implementation. The power-of-two boundary is checked separately.",
            body,
        ),
        PageBreak(),
        Paragraph("2. The fixed-modulus dimension boundary", heading),
        Paragraph(
            "Let <font name='Courier'>N=2^ell</font>. For constant <font name='Courier'>m,t</font>, using Theorem 3.1 as an ordinary clean-space algorithm would require <font name='Courier'>d=O(ell)</font> to reach <font name='Courier'>O(log n)</font> total space. Bhowmick-Dvir-Lovett Theorem 2, together with the now-proved bounded-torsion PFR theorem of Gowers-Green-Manners-Tao, gives for fixed <font name='Courier'>m</font>",
            body,
        ),
        Paragraph("MV(m,d) &lt;= exp(c_m d/log d).", formula),
        Paragraph(
            "Thus <font name='Courier'>N=2^ell</font> forces <font name='Courier'>d=Omega_m(ell log ell)</font>. This is consistent with the <font name='Courier'>d=O(ell log ell)</font> target in HPR Remark 3.4, which would give polynomial-time <font name='Courier'>O(log n log log n)</font> space, but it rules out reaching L through a fixed-modulus explicit family alone.",
            body,
        ),
        Paragraph(
            "<b>Question.</b> My reading is that the remaining L route must use the vectors through a succinct representation rather than materializing all <font name='Courier'>d</font> coordinates. Is that also how you view the boundary?",
            body,
        ),
        Paragraph("3. A characteristic-local selector variant", heading),
        Paragraph(
            "There is a small simplification if the update is projected to one CRT component <font name='Courier'>p_k</font>. At that component, use the shifted product itself instead of selecting one of its monomials. Its signed mixed difference is <font name='Courier'>xy-(x+1)y-x(y+1)+(x+1)(y+1)=1</font>.",
            body,
        ),
        Paragraph(
            "For <font name='Courier'>i!=k</font>, keep HPR's <font name='Courier'>(alpha_i,beta_i)</font> filters, normalize by their coefficients, and multiply by the CRT idempotent for <font name='Courier'>p_k</font>. Lemma 3.8 can be run modulo <font name='Courier'>m/p_k</font>, so the baseline inner products stored across its four oracle calls also omit the <font name='Courier'>p_k</font> component. The resulting one-level selector is one at the target and zero elsewhere while suspending <font name='Courier'>O(t+log(m/p_k))</font> bits rather than <font name='Courier'>O(log m)</font>.",
            body,
        ),
        Paragraph(
            "This identity is confined to its output characteristic. A fixed weight <font name='Courier'>phi:F_p-&gt;F_q</font> cannot provide the same mixed difference for distinct primes: summing over <font name='Courier'>x in F_p</font> makes the left side telescope to zero while the right side is <font name='Courier'>p!=0</font> in <font name='Courier'>F_q</font>. The positive identity passes 7,200 direct and 7,200 independently factorized checks over the <font name='Courier'>Z_15</font> companion family.",
            body,
        ),
        Paragraph(
            "I do not see a way to compose the projected updates recursively without recovering the complementary state, so this is only a one-level observation. <b>Question.</b> Could it still be useful for an asymmetric modulus, or is this limitation already captured by the CIR formulation?",
            body,
        ),
        Paragraph("References", heading),
        Paragraph("[1] A. Henzinger, E. Pyne, S. Ragavan. <i>Catalytic Tree Evaluation From Matching Vectors.</i> ECCC TR26-022, 2026.", refs),
        Paragraph("[2] A. Bhowmick, Z. Dvir, S. Lovett. <i>New Bounds for Matching Vector Families.</i> SIAM J. Comput. 43(5), 2014.", refs),
        Paragraph("[3] W. T. Gowers, B. Green, F. Manners, T. Tao. <i>Marton's Conjecture in Abelian Groups with Bounded Torsion.</i> 2024/2025.", refs),
        Paragraph("[4] H. Buhrman, R. Cleve, M. Koucky, B. Loff, F. Speelman. <i>Computing with a Full Memory: Catalytic Space.</i> STOC 2014.", refs),
        Paragraph("[5] W. Hesse, E. Allender, D. A. M. Barrington. <i>Uniform Constant-Depth Threshold Circuits for Division and Iterated Multiplication.</i> JCSS 65(4), 2002.", refs),
        Spacer(1, 3),
        Paragraph(
            "Implementation and full derivations are available if useful. Withdrawn ECCC TR26-044 is not used.",
            ParagraphStyle("Closing", parent=body, fontSize=7.8, leading=9.4, textColor=colors.HexColor("#596273")),
        ),
    ]
    doc.build(story)
    print(OUTPUT)


if __name__ == "__main__":
    build()
