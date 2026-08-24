/*
 * Swing pieces shared by the plugins that build their own windows.
 *
 * There are only two of those - smFRET Trace Histograms and smFRET Batch Analysis - and what they
 * share is one thing: a panel you can drop files onto. It is here rather than copied into both
 * because the rule about *where* to install it is not obvious, and a comment that only exists in
 * one of two copies is a comment that will be wrong in the other one before long.
 */

import javax.swing.TransferHandler;

import java.awt.datatransfer.DataFlavor;
import java.io.File;
import java.util.List;
import java.util.function.BooleanSupplier;
import java.util.function.Consumer;


final class smFRETSwing {

    private smFRETSwing() {
    }

    /**
     * A handler that accepts dropped files, for a panel that always wants them.
     */
    static TransferHandler fileDropHandler(Consumer<List<File>> onDrop) {
        return fileDropHandler(() -> true, onDrop);
    }

    /**
     * A handler that accepts dropped files while {@code accepting} says so.
     *
     * <p><b>Install it on the outermost container, once.</b> Not on the rows, and not again when
     * they are rebuilt. AWT picks the drop target by walking the component tree at the cursor for
     * the *deepest* component that has a DropTarget and falling back up the chain when a
     * descendant has none - see {@code Container.getMouseEventTarget} and
     * {@code DropTargetEventTargetFilter}. So a handler on a scroll pane covers its viewport, its
     * view, every row in the view and the blank space beside them, and rows built after the fact
     * are covered because nothing about them was ever registered.
     *
     * <p>The one thing that does intercept is a descendant with a DropTarget of its own, which in
     * practice means a text component - a JTextField accepts dropped text, so a file dropped on
     * one goes to it rather than here. Nothing in either of these panels is a text component.
     *
     * <p>{@code accepting} is asked on every drag over the panel, not once, so it can track state
     * that changes while the window is open - a batch that is running has nothing to add files to.
     *
     * @param accepting whether a drop should be taken right now
     * @param onDrop given the dropped files, on the event thread. Responsible for reporting its
     *        own failures: anything thrown out of it is taken as "the drop did not work" and
     *        nothing more, which is all this can honestly say about it.
     */
    static TransferHandler fileDropHandler(BooleanSupplier accepting, Consumer<List<File>> onDrop) {
        return new TransferHandler() {

            @Override
            public boolean canImport(TransferSupport support) {

                // isDrop() excludes a paste, which this would otherwise claim and then never be
                // sent - nothing here binds a paste key.
                return support.isDrop()
                        && support.isDataFlavorSupported(DataFlavor.javaFileListFlavor)
                        && accepting.getAsBoolean();
            }

            @Override
            @SuppressWarnings("unchecked")
            public boolean importData(TransferSupport support) {
                if (!canImport(support)) {
                    return false;
                }
                try {
                    onDrop.accept((List<File>) support.getTransferable()
                            .getTransferData(DataFlavor.javaFileListFlavor));
                    return true;
                } catch (Exception e) {
                    return false;
                }
            }
        };
    }
}
